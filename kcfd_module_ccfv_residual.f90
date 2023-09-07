module module_ccfv_residual
    implicit none
    
    public :: compute_residual
contains
    subroutine compute_residual
        
        use module_common_data     , only : p2, zero, half, one, two, four, three

        use module_common_data     , only : x, y, z, bc_type, nb
       
        use module_ccfv_data_grid  , only : nfaces, face, ncells, cell, &
                    bound, face_nrml, face_nrml_mag, face_centroid
       
        use module_ccfv_data_soln  , only : res, u, gradw, wsn, gradT, mu, w, gamma, u2w, Temp
        use module_input_parameter , only : Pr, Freestream_Temp ! Prandtl number
        use module_flux            , only : interface_flux, compute_viscosity, compute_tau, viscous_flux, viscous_bflux, &
                                            viscous_alpha
        use module_bc_states       , only : get_right_state, get_viscous_right_state
        use module_ccfv_gradient   , only : compute_gradient, compute_temperature_gradient
        use module_input_parameter , only : accuracy_order, use_limiter, navier_stokes
        use module_ccfv_limiter    , only : compute_limiter, phi
        use module_gewp            , only : gewp_solve

        implicit none

        real(p2)                    :: xm, ym, zm
        integer                     :: i, os
        integer                     :: c1, c2

        real(p2), dimension(5)      :: u1, u2
        real(p2), dimension(3,5)    :: gradw1, gradw2, gradwb
        real(p2), dimension(3)      :: unit_face_normal, bface_centroid

        real(p2), dimension(5)      :: num_flux
        integer                     :: j, ib, v1, v2, v3, ix, iu
        real(p2)                    :: xc1,xc2,yc1,yc2,zc1,zc2
        real(p2), dimension(5)      :: ub
        real(p2)                    :: wave_speed
        real(p2)                    :: phi1, phi2

        ! face_gradw = [dudx dudy dudz
        !               dvdx dvdy dvdz
        !               dwdx dwdy dwdz]
        real(p2), dimension(3,3)    :: face_gradw, gradQ, tau
        real(p2), dimension(3)      :: face_gradT, ds, dsds2, delU, w_face, grad_Tb
        ! face_gradT = [dTdx dTdy dTdz]
        real(p2)                    :: mu_face, heat_conductivity, rho_face
        real(p2), dimension(5)      :: visc_flux, wb, w1
        real(p2)                    :: theta_1, theta_2, theta_3

        real(p2), dimension(3)      :: d_Cb, V_parallel, wall_shear, cg_center, ejk
        real(p2)                    :: d_perp, magds, mag_ejk, Tb, a2
        real(p2), dimension(3,3)    :: del, del_inv
        integer                     :: bn1, bn2

        ! Debug
        ! real(p2), dimension(143)    :: Temp_local


        !--------------------------------------------------------------------------------
        ! Initialize the residuals and wsn = the sum of (max_wave_speed)*(face length)).
        ! Note: wsn is required to define a time step.
        cell_loop1 :  do i = 1, ncells

            res(:,i) = zero
            wsn(i)   = zero
     
        end do cell_loop1

        !--------------------------------------------------------------------------------
        ! Compute gradients at cells.

        !For now, let's set gradients to be zero in all cells...
        !Later we'll implement a subroutine that computes gradients.
        cell_loop2 : do i = 1, ncells
            gradw(:,:,i) = zero
        end do cell_loop2
        
        if (accuracy_order == 2 .or. navier_stokes) call compute_gradient
        if (navier_stokes) then
            call compute_temperature_gradient
            call compute_viscosity
        end if
        ! Temp_local = temp
        !local_gradw = gradw
        if (use_limiter)  call compute_limiter
        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        ! Residual computation: interior faces

        !--------------------------------------------------------------------------------
        ! Flux computation across internal faces (to be accumulated in res(:))
        !
        !          v2=Left(2)
        !        o---o---------o       face(j,:) = [i,k,v2,v1]
        !       .    .          .
        !      .     .           .
        !     .      .normal      .
        !    .  Left .--->  Right  .
        !   .   c1   .       c2     .
        !  .         .               .
        ! o----------o----------------o
        !          v1=Right(1)
        !
        !
        ! 1. Extrapolate the solutions to the face-midpoint from centroids 1 and 2.
        ! 2. Compute the numerical flux.
        ! 3. Add it to the residual for 1, and subtract it from the residual for 2.
        !
        !--------------------------------------------------------------------------------
        loop_faces : do i = 1,nfaces
            c1 = face(1,i)
            c2 = face(2,i)
            u1 = u(1:5,c1)
            u2 = u(1:5,c2)
            gradw1 = gradw(1:3,1:5,c1)
            gradw2 = gradw(1:3,1:5,c2)

            unit_face_normal = face_nrml(1:3,i)

            !use_limiter = .false. ! not implemented yet...
            if (use_limiter) then
                phi1 = phi(c1)
                phi2 = phi(c2)
            else 
                phi1 = one
                phi2 = one
            end if

            call interface_flux(          u1,       u2   , & !<- Left/right states
                                      gradw1,      gradw2, & !<- Left/right gradients
                                         unit_face_normal, & !<- unit face normal
                    cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                    cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, & !<- Right cell centroid
                                       face_centroid(1,i), &
                                       face_centroid(2,i), &
                                       face_centroid(3,i), & !<- face midpoint
                                        phi1,        phi2, & !<- Limiter functions
                                     num_flux, wave_speed  ) !<- Output)
            res(:,c1) = res(:,c1) + num_flux * face_nrml_mag(i)
            wsn(c1)   = wsn(c1) + wave_speed*face_nrml_mag(i)
            
            res(:,c2) = res(:,c2) - num_flux * face_nrml_mag(i)
            wsn(c2)   = wsn(c2) + wave_speed*face_nrml_mag(i)
            
            if (.not. navier_stokes) cycle

            ! ---------------------------------------------------------------
            ! Viscous Flux Terms
            ! ---------------------------------------------------------------

            ! ! Calculate the face gradients
            ! xc1 = cell(c1)%xc
            ! yc1 = cell(c1)%yc
            ! zc1 = cell(c1)%zc
            ! xc2 = cell(c2)%xc
            ! yc2 = cell(c2)%yc
            ! zc2 = cell(c2)%zc
            ! ds = (/xc2-xc1, yc2-yc1, zc2-zc1/) ! vector pointing from center of cell 1 to cell 2
            ! dsds2 = ds/(ds(1)**2 + ds(2)**2 + ds(3)**2) ! ds/ds^2
            ! do iu = 1,3
            !     delU = half * (gradw(:,iu+1,c1) + gradw(:,iu+1,c2))
            !     face_gradw(iu,:) = delU + ((w(iu+1,c2) - w(iu+1,c1)) - dot_product(delU,ds) ) * dsds2
            !     ! delU = delU_bar + [dU - dot(delU_bar,ds)]*ds/ds^2
            !     ! delU_bar is the arithmetic mean of the left and right gradients.  In order to prevent checkerboarding on certain 
            !     ! grids the gradient along the vector ds is replaced with a central difference.
            ! end do
            
            
            ! delU = half * (gradT(:,c1) + gradT(:,c2))
            ! face_gradT(:) = delU + ((Temp(c2) - Temp(c1)) - dot_product(delU,ds) ) * dsds2
            ! mu_face = half * (mu(c1) + mu(c2)) ! mu = M_inf*mu_ND/Re_inf
            ! tau = compute_tau(face_gradw,mu_face)
            ! heat_conductivity = mu_face/((gamma - one) * Pr)
            ! w_face = half * (w(2:4,c1) + w(2:4,c2))

            ! theta_1 = w_face(1)*tau(1,1) + w_face(2)*tau(2,1) + w_face(3)*tau(3,1) + heat_conductivity * face_gradT(1)
            ! theta_2 = w_face(1)*tau(1,2) + w_face(2)*tau(2,2) + w_face(3)*tau(3,2) + heat_conductivity * face_gradT(2)
            ! theta_3 = w_face(1)*tau(1,3) + w_face(2)*tau(2,3) + w_face(3)*tau(3,3) + heat_conductivity * face_gradT(3)


            ! visc_flux(1) = zero
            ! visc_flux(2) = face_nrml(1,i)*tau(1,1) + face_nrml(2,i)*tau(1,2) + face_nrml(3,i)*tau(1,3)
            ! visc_flux(3) = face_nrml(1,i)*tau(2,1) + face_nrml(2,i)*tau(2,2) + face_nrml(3,i)*tau(2,3)
            ! visc_flux(4) = face_nrml(1,i)*tau(3,1) + face_nrml(2,i)*tau(3,2) + face_nrml(3,i)*tau(3,3)
            ! visc_flux(5) = face_nrml(1,i)*theta_1  + face_nrml(2,i)*theta_2  + face_nrml(3,i)*theta_3

            ! write(*,*) visc_flux
            call viscous_flux(w(:,c1),w(:,c2), gradw1,gradw2, Temp(c1),Temp(c2), gradT(:,c1),gradT(:,c2), mu(c1),mu(c2), &
                                                 unit_face_normal, & !<- unit face normal
                            cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                            cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, & !<- Right cell centroid
                                                        visc_flux  )

            ! write(*,*) visc_flux
            ds = (/cell(c2)%xc - cell(c1)%xc,cell(c2)%yc-cell(c1)%yc,cell(c2)%zc-cell(c1)%zc/)
            magds = sqrt(ds(1)**2 + ds(2)**2 + ds(3)**2)
            ds = ds/magds
            ! call viscous_alpha(u1,u2,gradw1,gradw2,unit_face_normal,ds,magds,visc_flux)
            
            ! write(*,*) i,c1,c2,visc_flux
            ! ! if (isnan(visc_flux(3))) then
            ! !     write (*,*) "nan"
            ! ! end if
            ! write(*,*)
            ! Add the viscous flux
            res(:,c1) = res(:,c1) + visc_flux * face_nrml_mag(i)
            res(:,c2) = res(:,c2) - visc_flux * face_nrml_mag(i)
            
            ! Add to the spectral radius of the cell (CFD Principles and Applications Eq. 6.21)
            rho_face = half * (u(1,c1) + u(1,c2))
            wsn(c1) = wsn(c1) + (four/(cell(c1)%vol**2)) * max(four/(three*rho_face),gamma/rho_face) / (magds**0) & 
                        * two * (mu(c1)/Pr) * face_nrml_mag(i)
            wsn(c2) = wsn(c2) + (four/(cell(c2)%vol**2)) * max(four/(three*rho_face),gamma/rho_face) / (magds**0) & 
                        * two * (mu(c2)/Pr) * face_nrml_mag(i)


        end do loop_faces
        ! End of Residual computations for internal faces

        ! Now boundary faces
        boundary_loop : do ib = 1,nb
            bface_loop : do j = 1,bound(ib)%nbfaces
                bface_centroid = bound(ib)%bface_center(:,j)
                if (use_limiter) then
                    phi1 = phi(c1)
                    phi2 = one
                else 
                    phi1 = one
                    phi2 = one
                end if

                c1 = bound(ib)%bcell(j)
                u1 = u(1:5,c1)
                gradw1 = gradw(1:3,1:5,c1)
                unit_face_normal = bound(ib)%bface_nrml(:,j)

                ! Get the right hand state (weak BC!)
                call get_right_state(bface_centroid(1), bface_centroid(2), bface_centroid(3), &
                                     u1, unit_face_normal, bc_type(ib), ub)
                gradw2 = zero ! won't matter since boundary cell center will be at face center

                call interface_flux(          u1,      ub, & !<- Left/right states
                                      gradw1,      gradw2, & !<- Left/right gradients
                                         unit_face_normal, & !<- unit face normal
                    cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                                        bface_centroid(1), &
                                        bface_centroid(2), &
                                        bface_centroid(3), & !<- Right cell centroid
                                        bface_centroid(1), &
                                        bface_centroid(2), &
                                        bface_centroid(3), & !<- boundary ghost cell "center"
                                        phi1,        phi2, & !<- Limiter functions
                                     num_flux, wave_speed  ) !<- Output)
                res(:,c1) = res(:,c1) + num_flux * bound(ib)%bface_nrml_mag(j)
                wsn(c1)   = wsn(c1) + wave_speed * bound(ib)%bface_nrml_mag(j)         
                !local_res = res

                if (.not. navier_stokes) then 
                    cycle
                    ! Add wall shear forces
                else
                    ! This computation is based off of section 15.6.1.1 in "The Finite Volume Method in Computational 
                    ! Fluid Dynamics" (url: https://link.springer.com/book/10.1007/978-3-319-16874-6)

                    ! wb = u2w(ub)

                    ! mu_face = mu(c1) ! just set it to the cell value.  This should be sufficient.  We can be more clever later if we
                    ! ! want...
                    
                    ! ! Vector from cell center to boundary face center: 
                    ! d_Cb = -(/cell(c1)%xc - bface_centroid(1), cell(c1)%yc - bface_centroid(2), &
                    !          cell(c1)%zc - bface_centroid(3)/)
                    ! ! Component of d_Cb perpendicular to boundary face:
                    ! d_perp = dot_product(d_Cb,unit_face_normal)

                    ! ! Component of velocity at cell center that is parallel to boundary face
                    ! ! V_par = V - (V * n)n
                    ! ! V_parallel = w(2:4,c1) - dot_product(w(2:4,c1),unit_face_normal) * unit_face_normal

                    ! ! Calculate the wall shear force:
                    ! wall_shear = -(mu_face*bound(ib)%bface_nrml_mag(j)/d_perp) * &
                    !             ((w(2:4,c1) - wb(2:4)) - (dot_product(w(2:4,c1) - wb(2:4),unit_face_normal)*unit_face_normal))
                    ! ! I think I have the sign right...
                    ! res(2:4,c1) = res(2:4,c1) + wall_shear * bound(ib)%bface_nrml_mag(j)

                    ! ! Calculating wall temperature gradient
                    
                    ! ! Calculate two vectors from boundary face center to boundary nodes.
                    ! bn1 = bound(ib)%bfaces(2,j)
                    ! bn2 = bound(ib)%bfaces(3,j)
                    ! del(1,:) = d_Cb
                    ! del(2,:) = (/x(bn1), y(bn1), z(bn1)/) - bface_centroid
                    ! del(3,:) = (/x(bn2), y(bn2), z(bn2)/) - bface_centroid

                    ! call gewp_solve(del, 3, del_inv, os)

                    ! face_gradT = matmul(del_inv,(/one - Temp(c1), zero, zero/))
                    ! ! write(*,*) gradT(:,c1), ',', Temp(c1)
                    ! heat_conductivity = mu_face/((gamma - one) * Pr) ! T_wall = T_inf = 1 (T is normalized with T_inf)

                    ! res(5,c1) = res(5,c1) - heat_conductivity *  dot_product(face_gradT,unit_face_normal) &
                    !             * bound(ib)%bface_nrml_mag(j)

                    ! call viscous_bflux(w(:,c1),wb,Temp(c1),one,mu(c1), &
                    !           cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                    !                               bface_centroid(1), &
                    !                               bface_centroid(2), &
                    !                               bface_centroid(3), & 
                    !                           bound(ib)%bfaces(2,j), &
                    !                           bound(ib)%bfaces(3,j), &
                    !                     bound(ib)%bface_nrml_mag(j), &
                    !                       bound(ib)%bface_nrml(:,j), &
                    !                                       visc_flux  )

                    !* Building a ghost cell
                    !*           njk
                    !*  Face normal ^   o Ghost cell center
                    !*              |  .
                    !*              | . ejk (edge vector from boundary cell to ghost cell)
                    !*              |.
                    !*       -------x-------- Face
                    !*             .
                    !*            . d_Cb (vector from boundary cell center to boundary face center = 0.5*ejk)
                    !*           .
                    !*          o Boundary cell center

                    d_Cb = (/bface_centroid(1) - cell(c1)%xc,bface_centroid(2) - cell(c1)%yc,bface_centroid(3) - cell(c1)%zc/)
                    ejk = 2*d_Cb
                    mag_ejk = sqrt( ejk(1)**2 + ejk(2)**2 + ejk(3)**2 )
                    ejk = ejk/mag_ejk
                    cg_center = bface_centroid + d_Cb ! ghost cell center
                    
                    
                    call get_viscous_right_state(bface_centroid(1),bface_centroid(2),bface_centroid(3), &
                                                 u1,gradw1,unit_face_normal,bc_type(ib),ub,gradwb)

                    ! call viscous_alpha(u1,ub,gradw1,gradwb,bound(ib)%bface_nrml(:,j),ejk,mag_ejk,visc_flux)
                    ! write(*,*) visc_flux
                    w1 = w(:,c1)
                    wb = u2w(ub)
                    Tb = gamma*wb(5)/wb(1)
                    a2 = gamma*half*(wb(5))/wb(1)
                    grad_Tb = ( gamma*gradwb(:,5) - a2*gradwb(:,1)) /wb(1)

                    call viscous_flux(w1,wb, gradw1,gradwb, Temp(c1),Tb, gradT(:,c1),grad_Tb, mu(c1),mu(c1), &
                                                 unit_face_normal, & !<- unit face normal
                          cell(c1)%xc,  cell(c1)%yc,  cell(c1)%zc, & !<- Left  cell centroid
                         cg_center(1), cg_center(2), cg_center(3), & !<- Right cell centroid
                                                        visc_flux  )

                    res(:,c1) = res(:,c1) + visc_flux * bound(ib)%bface_nrml_mag(j)
                    if (any(isnan(res(:,c1)))) then 
                        write (*,*) "nan value present - press [Enter] to continue"
                        read(unit=*,fmt=*)
                    end if
                    wsn(c1) = wsn(c1) + (four/(cell(c1)%vol**2)) * max(four/(three*u(1,c1)),gamma/u(1,c1))/(mag_ejk**0) &
                                * two * (mu(c1)/Pr) * bound(ib)%bface_nrml_mag(j)
                end if
            end do bface_loop
        end do boundary_loop

    end subroutine compute_residual
end module module_ccfv_residual