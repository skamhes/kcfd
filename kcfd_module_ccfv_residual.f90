module module_ccfv_residual
    implicit none
    
    public :: compute_residual
contains
    subroutine compute_residual
        
        use module_common_data     , only : p2, zero, half, one, two, four, three, third

        use module_common_data     , only : x, y, z, bc_type, nb
       
        use module_ccfv_data_grid  , only : nfaces, face, ncells, cell, &
                    bound, face_nrml, face_nrml_mag, face_centroid
       
        use module_ccfv_data_soln  , only : res, u, gradw, wsn, gradT, mu, w, gamma, u2w, Temp, uR2, gradq, q, u2q, q2u
        use module_input_parameter , only : Pr, Freestream_Temp, visc_flux_method ! Prandtl number
        use module_flux            , only : interface_flux, compute_viscosity, compute_tau, viscous_flux, viscous_bflux, &
                                            viscous_alpha
        use module_bc_states       , only : get_right_state, get_viscous_right_state
        use module_ccfv_gradient   , only : compute_gradient, compute_temperature_gradient
        use module_input_parameter , only : accuracy_order, use_limiter, navier_stokes, low_mach_correction
        use module_ccfv_limiter    , only : compute_limiter, phi
        use module_gewp            , only : gewp_solve

        implicit none

        real(p2)                    :: xm, ym, zm
        integer                     :: i, os
        integer                     :: c1, c2

        real(p2), dimension(5)      :: u1, u2, q1, q2
        real(p2), dimension(3,5)    :: gradw1, gradw2, gradwb, gradq1, gradq2, gradqb
        real(p2), dimension(3)      :: unit_face_normal, bface_centroid

        real(p2), dimension(5)      :: num_flux
        integer                     :: j, ib, v1, v2, v3, ix, iu
        real(p2)                    :: xc1,xc2,yc1,yc2,zc1,zc2
        real(p2), dimension(5)      :: ub, qb
        real(p2)                    :: uR2b
        real(p2)                    :: wave_speed
        real(p2)                    :: phi1, phi2

        ! face_gradw = [dudx dudy dudz
        !               dvdx dvdy dvdz
        !               dwdx dwdy dwdz]
        real(p2), dimension(3,3)    :: face_gradw, tau
        real(p2), dimension(3)      :: face_gradT, ds, dsds2, delU, w_face, grad_Tb
        ! face_gradT = [dTdx dTdy dTdz]
        real(p2)                    :: mu_face, heat_conductivity, rho_face
        real(p2), dimension(5)      :: visc_flux, wb, w1
        real(p2)                    :: theta_1, theta_2, theta_3

        real(p2), dimension(3)      :: d_Cb, V_parallel, wall_shear, cg_center, ejk
        real(p2)                    :: d_perp, magds, mag_ejk, Tb, a2, absU, a
        real(p2), dimension(3,3)    :: del, del_inv
        integer                     :: bn1, bn2

        real(p2) :: epsilon = 1e-05_p2
        real(p2) :: local_vel

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
        if (low_mach_correction) then
            cell_loop2_q : do i = 1, ncells
                gradq(:,:,i) = zero
            end do cell_loop2_q
        else
            cell_loop2_w : do i = 1, ncells
                gradw(:,:,i) = zero
            end do cell_loop2_w
        end if
        
        if (low_mach_correction) then
            do i = 1,ncells
                absU = sqrt(q(2,i)**2 + q(3,i)**2 + q(4,i)**2)
                a = q(5,i)
                if (absU < epsilon * a) then
                    uR2(i) = (epsilon * a) ** 2
                elseif (absU < a) then
                    uR2(i) = absU ** 2
                else
                    uR2(i) = a**2
                end if
            end do

            if (navier_stokes) then
                ! loop the faces.  I'm sure this can be done more efficiently but right now im just working on getting it working...
                do i = 1,nfaces
                    c1 = face(1,i)
                    c2 = face(2,i)
                    ejk = (/ cell(c2)%xc - cell(c1)%xc, cell(c2)%yc - cell(c1)%yc, cell(c2)%zc - cell(c1)%zc /)
                    mag_ejk = sqrt(ejk(1)**2 + ejk(2)**2 + ejk(3)**2)
                
                    if ( sqrt(uR2(c1)) < (mu(c1)/(w(1,c1)*mag_ejk)) ) then
                        uR2(c1) = (mu(c2)/(w(1,c2)*mag_ejk)) ** 2
                    end if
                    if ( sqrt(uR2(c2)) < (mu(c2)/(w(1,c2)*mag_ejk)) ) then
                        uR2(c2) = (mu(c2)/(w(1,c2)*mag_ejk)) ** 2
                    end if
                end do
            end if                
        end if

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
            if (low_mach_correction) then
                q1 = q(1:5, c1)
                q2 = q(1:5, c2)
                gradq1 = gradq(1:3,1:5,c1)
                gradq2 = gradq(1:3,1:5,c2)

            else
                u1 = u(1:5,c1)
                u2 = u(1:5,c2)
                gradw1 = gradw(1:3,1:5,c1)
                gradw2 = gradw(1:3,1:5,c2)
            end if

            unit_face_normal = face_nrml(1:3,i)

            !use_limiter = .false. ! not implemented yet...
            if (use_limiter) then
                phi1 = phi(c1)
                phi2 = phi(c2)
            else 
                phi1 = one
                phi2 = one
            end if
            if (low_mach_correction) then
                
                call interface_flux(          q1,       q2   , & !<- Left/right states
                                          gradq1,      gradq2, & !<- Left/right gradients
                                             unit_face_normal, & !<- unit face normal
                        cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                        cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, & !<- Right cell centroid
                                           face_centroid(1,i), &
                                           face_centroid(2,i), &
                                           face_centroid(3,i), & !<- face midpoint
                                            phi1,        phi2, & !<- Limiter functions
                                         num_flux, wave_speed, & !<- Output)
                                     uR2L=uR2(c1), uR2R=uR2(c2)) !<- scaling factor (optional)
            else
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
            end if
            
            res(:,c1) = res(:,c1) + num_flux * face_nrml_mag(i)
            wsn(c1)   = wsn(c1) + wave_speed*face_nrml_mag(i)
            
            res(:,c2) = res(:,c2) - num_flux * face_nrml_mag(i)
            wsn(c2)   = wsn(c2) + wave_speed*face_nrml_mag(i)
            
            if (.not. navier_stokes) cycle
            
            ds = (/cell(c2)%xc - cell(c1)%xc,cell(c2)%yc-cell(c1)%yc,cell(c2)%zc-cell(c1)%zc/)
            magds = sqrt(ds(1)**2 + ds(2)**2 + ds(3)**2)
            ds = ds/magds
            if (visc_flux_method == 'corrected') then ! Corrected gradient
                call viscous_flux(w(:,c1),w(:,c2), gradw1,gradw2, Temp(c1),Temp(c2), gradT(:,c1),gradT(:,c2), mu(c1),mu(c2), &
                                                    unit_face_normal, & !<- unit face normal
                                cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                                cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, & !<- Right cell centroid
                                                            visc_flux  )
            else if (visc_flux_method == 'alpha') then
                call viscous_alpha(u1,u2,gradw1,gradw2,unit_face_normal,ds,magds,visc_flux)
            else
                write(*,*) 'Unsupported viscous flux method. Stop!'
                stop
            end if

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
                
                
                unit_face_normal = bound(ib)%bface_nrml(:,j)

                
                gradw2 = zero ! won't matter since boundary cell center will be at face center

                if (low_mach_correction) then
                    q1 = q(:,c1)
                    u1 = q2u(q1)
                    ! Get the right hand state (weak BC!)
                    call get_right_state(bface_centroid(1), bface_centroid(2), bface_centroid(3), &
                                         u1, unit_face_normal, bc_type(ib), ub)
                    gradq1 = gradq(1:3,1:5,c1)
                    qb = u2q(ub)
                    uR2b = min(qb(5) , max(qb(2)**2 + qb(3)**2 + qb(4)**2, epsilon*qb(5) ))
                    call interface_flux(          q1,      qb, & !<- Left/right states
                                          gradq1,      gradw2, & !<- Left/right gradients
                                             unit_face_normal, & !<- unit face normal
                        cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                                            bface_centroid(1), &
                                            bface_centroid(2), &
                                            bface_centroid(3), & !<- Right cell centroid
                                            bface_centroid(1), &
                                            bface_centroid(2), &
                                            bface_centroid(3), & !<- boundary ghost cell "center"
                                            phi1,        phi2, & !<- Limiter functions
                                         num_flux, wave_speed, & !<- Output)
                                     uR2L=uR2(c1), uR2R=uR2b) !<- scaling factor (optional)
                else
                    u1 = u(1:5,c1)
                    ! Get the right hand state (weak BC!)
                    call get_right_state(bface_centroid(1), bface_centroid(2), bface_centroid(3), &
                                         u1, unit_face_normal, bc_type(ib), ub)
                    gradw1 = gradw(1:3,1:5,c1)
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
                end if
                res(:,c1) = res(:,c1) + num_flux * bound(ib)%bface_nrml_mag(j)
                wsn(c1)   = wsn(c1) + wave_speed * bound(ib)%bface_nrml_mag(j)         
                !local_res = res

                if (.not. navier_stokes) then 
                    cycle
                    ! Add wall shear forces
                else
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

                    if (visc_flux_method == 'alpha') then
                        call viscous_alpha(u1,ub,gradw1,gradwb,bound(ib)%bface_nrml(:,j),ejk,mag_ejk,visc_flux)
                    elseif (visc_flux_method == 'corrected') then
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
                    else
                        write(*,*) 'Unsupported viscous flux method. Stop!'
                        stop
                    end if

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