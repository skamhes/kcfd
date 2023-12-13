module module_ccfv_residual
    implicit none
    
    public :: compute_residual
contains
    subroutine compute_residual
        
        use module_common_data     , only : p2, zero, half, one, two, four, three, third

        use module_common_data     , only : x, y, z, bc_type, nb
       
        use module_ccfv_data_grid  , only : nfaces, face, ncells, cell, &
                    bound, face_nrml, face_nrml_mag, face_centroid
       
        use module_ccfv_data_soln  , only : res, u, gradw, wsn, gradT, mu, w, gamma, u2w, Temp, uR2, gradq, q, u2q, q2u, &
                                            gradq_w, gradw_w
        use module_input_parameter , only : Pr, Freestream_Temp, visc_flux_method ! Prandtl number
        use module_flux            , only : interface_flux, compute_viscosity, compute_tau, viscous_flux, viscous_bflux, &
                                            viscous_alpha, viscous_flux_primitive
        use module_bc_states       , only : get_right_state, get_viscous_right_state
        use module_ccfv_gradient   , only : compute_gradient, compute_temperature_gradient, compute_weighted_gradient
        use module_input_parameter , only : accuracy_order, use_limiter, navier_stokes, low_mach_correction, min_ref_vel, &
                                            eps_weiss_smith, pressure_dif_ws
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

        real(p2) :: dp
        ! real(p2) :: epsilon = 1e-05_p2
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
                if (navier_stokes) gradq_w(:,:,i) = zero
            end do cell_loop2_q
        else
            cell_loop2_w : do i = 1, ncells
                gradw(:,:,i) = zero
                if (navier_stokes) gradw_w(:,:,i) = zero
            end do cell_loop2_w
        end if
        
        if (low_mach_correction) then
            do i = 1,ncells
                uR2(i) = min(q(5,i),max(sqrt(q(2,i)**2 + q(3,i)**2 + q(4,i)**2), eps_weiss_smith*q(5,i),min_ref_vel))**2
                ! uR2(i) = min(q(5,i),max(sqrt(q(2,i)**2 + q(3,i)**2 + q(4,i)**2),min_ref_vel))**2
            end do
            do i = 1,nfaces
                c1 = face(1,i)
                c2 = face(2,i)
                if (navier_stokes) then
                    !stuff
                end if
                dp = abs(q(1,c2) - q(1,c1))
                uR2(c1) = min(q(5,c1),max(sqrt(uR2(c1)), pressure_dif_ws*sqrt(dp*q(5,c1)/(gamma*q(1,c1))) ))**2
                uR2(c2) = min(q(5,c2),max(sqrt(uR2(c2)), pressure_dif_ws*sqrt(dp*q(5,c2)/(gamma*q(1,c2))) ))**2
            end do
        end if

        if (accuracy_order == 2) call compute_gradient
        if (navier_stokes) call compute_weighted_gradient
        if (navier_stokes) then
            if (.not. low_mach_correction) call compute_temperature_gradient
            call compute_viscosity
        end if

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
                if (low_mach_correction) then
                    call viscous_flux_primitive(q1,q2,gradq_w(:,:,c1),gradq_w(:,:,c2),mu(c1),mu(c2), &
                                                     unit_face_normal, & !<- unit face normal
                                cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                                cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, & !<- Right cell centroid
                                                            visc_flux  )
                else
                    call viscous_flux(w(:,c1),w(:,c2), gradw_w(:,:,c1),gradw_w(:,:,c2), Temp(c1),Temp(c2), gradT(:,c1),gradT(:,c2) &
                                                          , mu(c1),mu(c2), &
                                                         unit_face_normal, & !<- unit face normal
                                    cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                                    cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, & !<- Right cell centroid
                                                                visc_flux  )
                end if
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
            if (low_mach_correction) then
                rho_face = half*gamma*((q1(1)/q1(5)) + (q2(1)/q2(5)))
            else
                rho_face = half * (u1(1) + u2(1))
            end if
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
                    uR2b = min(qb(5) , max(sqrt(qb(2)**2 + qb(3)**2 + qb(4)**2), eps_weiss_smith*qb(5),min_ref_vel ))**2
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
                                         num_flux, wave_speed  , & !<- Output)
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
                    if (low_mach_correction) then
                        call get_viscous_right_state(bface_centroid(1),bface_centroid(2),bface_centroid(3), &
                                                     u1,gradq_w(:,:,c1),unit_face_normal,bc_type(ib),ub,gradqb)
                        qb = u2q(ub)
                    else
                        call get_viscous_right_state(bface_centroid(1),bface_centroid(2),bface_centroid(3), &
                                                     u1,gradw_w(:,:,c1),unit_face_normal,bc_type(ib),ub,gradwb)
                    end if
                    if (visc_flux_method == 'alpha') then
                        call viscous_alpha(u1,ub,gradw1,gradwb,bound(ib)%bface_nrml(:,j),ejk,mag_ejk,visc_flux)
                    elseif (visc_flux_method == 'corrected') then
                        if (low_mach_correction) then
                            call viscous_flux_primitive(q1,qb,gradq_w(:,:,c1),gradqb,mu(c1),mu(c1), &
                                                                         unit_face_normal, & !<- unit face normal
                                                  cell(c1)%xc,  cell(c1)%yc,  cell(c1)%zc, & !<- Left  cell centroid
                                                 cg_center(1), cg_center(2), cg_center(3), & !<- Right cell centroid
                                                                                visc_flux  )
                        else
                            w1 = w(:,c1)
                            wb = u2w(ub)
                            Tb = gamma*wb(5)/wb(1)
                            a2 = gamma*half*(wb(5))/wb(1)
                            grad_Tb = ( gamma*gradwb(:,5) - a2*gradwb(:,1)) /wb(1)

                            call viscous_flux(w1,wb, gradw_w(:,:,c1),gradwb, Temp(c1),Tb, gradT(:,c1),grad_Tb, mu(c1),mu(c1), &
                                                        unit_face_normal, & !<- unit face normal
                                cell(c1)%xc,  cell(c1)%yc,  cell(c1)%zc, & !<- Left  cell centroid
                                cg_center(1), cg_center(2), cg_center(3), & !<- Right cell centroid
                                                                visc_flux  )
                        end if
                    else
                        write(*,*) 'Unsupported viscous flux method. Stop!'
                        stop
                    end if

                    res(:,c1) = res(:,c1) + visc_flux * bound(ib)%bface_nrml_mag(j)
                    if (any(isnan(res(:,c1)))) then 
                        write (*,*) "nan value present cell = ",c1," - press [Enter] to continue"
                        read(unit=*,fmt=*)
                    end if
                    if (low_mach_correction) then
                        rho_face = gamma*(q1(1)/q1(5)) ! not actually the face, just reusing variables...
                    else
                        rho_face = (u1(1) + u2(1))
                    end if

                    wsn(c1) = wsn(c1) + (four/(cell(c1)%vol**2)) * max(four/(three*rho_face),gamma/rho_face)/(mag_ejk**0) &
                                * two * (mu(c1)/Pr) * bound(ib)%bface_nrml_mag(j)
                end if

            end do bface_loop
        end do boundary_loop

    end subroutine compute_residual

    subroutine solution_limiting(res,q_ref,grad_pressure,vol,dtau)
        
        ! subroutine for limiting the CFL in areas where large change is expected to prevent divergence
        ! https://doi.org/10.1016/j.jcp.2009.03.040

        ! I would have put this in the steady solver module but circular dependencies and stuff...  Plus this uses the residual!

        use module_common_data , only : p2, half, third
        use module_input_parameter , only : CFL, M_inf
        use module_ccfv_data_soln , only : gamma, p_inf

        implicit none

        real(p2), dimension(5), intent(in   ) :: res ! estimated solution update Q^{n+1} (equal to cell residual)
        real(p2), dimension(5), intent(in   ) :: q_ref ! current solution Q^{n}
        real(p2), dimension(3), intent(in   ) :: grad_pressure ! gradient of pressure in the cell
        real(p2),               intent(in   ) :: vol ! cell volume (used for approximating pressure diff accross cell)
        real(p2),               intent(inout) :: dtau ! (potentially) limited timestep

        real(p2) :: p_ref, T_ref         ! reference values that will be used to see if we need to limit the solution
        real(p2), dimension(3) :: u_ref  ! reference values that will be used to see if we need to limit the solution
        real(p2) :: alpha = 0.1_p2       ! maximum allowable fractional change
        real(p2) :: alpha_u = 2._p2      ! maximum allowable change for local velocity
        real(p2) :: p_grad_mag           ! magnitude of the pressure gradient
        real(p2) :: u_mag                ! magnitude of cell vellocity
        real(p2) :: delta_p              ! measurement of variation in p across the cell
        real(p2) :: min_limit = 1e-09_p2 ! minimum limit to prevent ref values going to zero
        real(p2) :: rho                  ! local cell density (wer're using primitive variables after all)
        real(p2) :: CFL_mod              ! limited CFL
        real(p2), dimension(5 ):: q_est  ! estimated change in the solution after update (q_est = dTau * res)
        integer :: i


        ! First Calculate rho.  We'll need it.
        rho = q_ref(1)*gamma/q_ref(5)

        ! Approximate the variation in pressure across the cell.  This is done by taking the cell "length" and multiplying by the 
        ! pressure gradient.  For 3D cells vol^1/3 is an ok approximation of the cell length.  For now we're gonna try multiplying 
        ! by a half to be a little more conservative.
        p_grad_mag = sqrt(grad_pressure(1)**2 + grad_pressure(2)**2 + grad_pressure(3)**2)
        delta_p = half * p_grad_mag * vol ** (third)

        ! Caluculate magnitude of cell velocity
        u_mag = sqrt(q_ref(2)**2 + q_ref(3)**2 + q_ref(4)**2)

        ! First calculate p_ref
        p_ref = alpha * max(half*rho*u_mag**2 , delta_p , p_inf*min_limit)

        ! Next the velocity ref for each direction
        do i = 1,3
            u_ref(i) = max(alpha_u*abs(q_ref(i+1) ), alpha * q_ref(5) * delta_p / (gamma * q_ref(1)) , M_inf * min_limit )
        end do

        ! The reference temperature is simple
        T_ref = max(alpha * q_ref(5), min_limit)

        ! Calculate q_est = dTau * res and check for alpha
        q_est = dtau * res

        ! First set the initial CFL
        CFL_mod = CFL

        ! Check pressure update
        CFL_mod = min( CFL_mod,CFL * p_ref / abs(q_est(1)) )

        ! repeat for velocity
        do i = 1,3
            CFL_mod = min( CFL_mod, CFL * u_ref(i) / abs(q_est(i+1)) )
        end do

        ! Repeat for temperature
        CFL_mod = min(CFL_mod, CFL * T_ref / abs(q_est(5)) )

        ! if the CFL has been limited apply the limiting to dtau
        if (CFL_mod < CFL) then
            dtau = dtau * (CFL_mod / CFL)
        end if

    end subroutine solution_limiting
end module module_ccfv_residual