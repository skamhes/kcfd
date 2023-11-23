module module_numerical_jacobian


    implicit none

contains

    subroutine compute_numerical_jacobian

        use module_common_data, only : p2, my_eps, one, canonical_array, zero, two, bc_type, nb, four, half
        use module_ccfv_data_grid       , only : nfaces, face, ncells, cell, face_nrml, face_nrml_mag, jac, &
                                                 kth_nghbr_of_1, kth_nghbr_of_2, bound, face_centroid
        use module_ccfv_data_soln       , only : w, u, u2w, dtau, Temp, mu, gradw, gradT, gamma, q, q2u, gradq, u2q, wsn, uR2
        use module_bc_states            , only : get_right_state
        use module_flux                 , only : interface_flux, perturb_flux, get_viscosity, viscous_flux, viscous_bflux
        use module_gewp                 , only : gewp_solve
        use module_input_parameter      , only : navier_stokes, use_limiter, low_mach_correction
        use module_ccfv_limiter         , only : phi
        use module_ccfv_gradient        , only : perturbation_gradient, perturbation_gradient_boundary, &
                                                 perturb_temperature_gradient
        
        
        implicit none
        integer                     :: c1, c2, i, k, ib, idestat, j, ivar, ii

        real(p2), dimension(5)      :: u1, u2, u_perturb, w1,w2, w_perturb, q1, q2
        real(p2), dimension(3,5)    :: gradw1, gradw2, grad_perturb1, grad_perturb2
        real(p2), dimension(3)      :: unit_face_normal, bface_centroid

        real(p2), dimension(5)      :: num_flux, base_num_flux, visc_flux, base_visc_flux
        integer                     :: v1, v2, v3, ix, iu
        real(p2)                    :: xc1,xc2,yc1,yc2,zc1,zc2
        real(p2), dimension(5)      :: ub, wb, qb
        real(p2)                    :: wave_speed, face_mag
        real(p2)                    :: phi1, phi2

        real(p2)                    :: T_perturb, mu_perturb
        real(p2), dimension(3)      :: gradT_perturb1,gradT_perturb2

        real(p2)                    :: perturbation_value

        real(p2), dimension(5,5)    :: dFnduL, dFnduR

        real(p2), dimension(5,5)    :: preconditioner
        real(p2), dimension(5,5)    :: duLdqL, duRdqR
        real(p2)                    :: theta
        real(p2)                    :: rho_p, rho_T, rho
        real(p2)                    :: H, alpha, beta, lambda, absu, uR2b
        real(p2) :: epsilon = 1e-05_p2

        perturbation_value = four * sqrt(my_eps) ! empirical decision
        


        ! Initialize jacobian terms
        do i = 1,ncells
            jac(i)%diag = zero
            jac(i)%off_diag = zero
            jac(i)%diag_inv = zero
        end do

        loop_faces : do i = 1,nfaces
            c1 = face(1,i)
            c2 = face(2,i)
            if (low_mach_correction) then
                q1 = q(:,c1)
                q2 = q(:,c2)
                u1 = q2u(q1)
                u2 = q2u(q2)
                gradw1 = gradq(:,:,c1)
                gradw2 = gradq(:,:,c2)
            else
                u1 = u(1:5,c1)
                u2 = u(1:5,c2)
                gradw1 = gradw(1:3,1:5,c1)
                gradw2 = gradw(1:3,1:5,c2)
            endif


            face_mag       = face_nrml_mag(i)
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
                                          gradw1,      gradw2, & !<- Left/right gradients
                                             unit_face_normal, & !<- unit face normal
                        cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                        cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, & !<- Right cell centroid
                                           face_centroid(1,i), &
                                           face_centroid(2,i), &
                                           face_centroid(3,i), & !<- face midpoint
                                            phi1,        phi2, & !<- Limiter functions
                                    base_num_flux, wave_speed, &!<- Output)
                                    uR2L=uR2(c1), uR2R=uR2(c2) ) 
            else ! Regular
                call interface_flux(          u1,       u2   , & !<- Left/right states
                                          gradw1,      gradw2, & !<- Left/right gradients
                                             unit_face_normal, & !<- unit face normal
                        cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                        cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, & !<- Right cell centroid
                                           face_centroid(1,i), &
                                           face_centroid(2,i), &
                                           face_centroid(3,i), & !<- face midpoint
                                            phi1,        phi2, & !<- Limiter functions
                                     base_num_flux, wave_speed  ) !<- Output)
            endif

            loop_sides : do ii = 1,2
                loop_var_perturbs : do ivar = 1,5
                    if (low_mach_correction) then
                        call perturb_flux(   q1,   q2, & !<- Left/right states
                                    gradw1,    gradw2, & !<- Left/right gradients
                                     unit_face_normal, & !<- unit face normal
                cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, & !<- Right cell centroid
                                   face_centroid(1,i), &
                                   face_centroid(2,i), &
                                   face_centroid(3,i), & !<- face midpoint
                                    phi1,        phi2, & !<- Limiter functions
                         ivar, perturbation_value, ii, & !<- Var to be perturbed, size of pert and side of face
                                 num_flux, wave_speed, & !<- Output)
                            uR2L=uR2(c1), uR2R=uR2(c2) )
                    else
                        call perturb_flux(   u1,   u2, & !<- Left/right states
                                    gradw1,    gradw2, & !<- Left/right gradients
                                     unit_face_normal, & !<- unit face normal
                cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, & !<- Right cell centroid
                                   face_centroid(1,i), &
                                   face_centroid(2,i), &
                                   face_centroid(3,i), & !<- face midpoint
                                    phi1,        phi2, & !<- Limiter functions
                         ivar, perturbation_value, ii, & !<- Var to be perturbed, size of pert and side of face
                                  num_flux, wave_speed  ) !<- Output)
                    end if
                    if (ii == 1) then ! left side
                        ! u_perturb = u1
                        ! u_perturb(ivar) = u1(ivar) + perturbation_value
                        ! call perturbation_gradient(gradw1,gradw2, ivar, c1,c2, u_perturb, grad_perturb1,grad_perturb2)
                        !
                        ! Now 
                        dFnduL(:,ivar) = (num_flux - base_num_flux) / perturbation_value
                    else
                        ! u_perturb = u2 
                        ! u_perturb(ivar) = u2(ivar) + perturbation_value
                        ! call perturbation_gradient(gradw2, gradw1, ivar, c2,c1, u_perturb, grad_perturb2,grad_perturb1)
                        !
                        ! Now 
                        dFnduR(:,ivar) = (num_flux - base_num_flux) / perturbation_value
                    end if
                end do loop_var_perturbs
            end do loop_sides
          
            ! Add to diagonal term of C1
            jac(c1)%diag            = jac(c1)%diag            + dFnduL * face_mag
            ! get neighbor index k for cell c1
            k = kth_nghbr_of_1(i)
            ! add to off diagonal neighbor k for cell c1
            jac(c1)%off_diag(:,:,k) = jac(c1)%off_diag(:,:,k) + dFnduR * face_mag

            ! Subtract terms from c2
            jac(c2)%diag            = jac(c2)%diag            - dFnduR * face_mag
            k = kth_nghbr_of_2(i)
            jac(c2)%off_diag(:,:,k) = jac(c2)%off_diag(:,:,k) - dFnduL * face_mag

            !-------------------------------
            ! Viscous Flux
            !-------------------------------

            ! note: this can be combined with the loops above but for now were gonna seperate the two for readability...
            w1 = u2w(u1)
            w2 = u2w(u2)

            if (.not. navier_stokes) cycle loop_faces

            call viscous_flux(w(:,c1),w(:,c2), gradw1,gradw2, Temp(c1),Temp(c2), gradT(:,c1),gradT(:,c2), mu(c1),mu(c2), &
                                                 unit_face_normal, & !<- unit face normal
                            cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                            cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, & !<- Right cell centroid
                                                   base_visc_flux  )

            loop_sides_v : do ii = 1,2
                loop_var_perturbs_v : do ivar = 1,5
                    if (ii == 1) then ! left side
                        ! Perturb U
                        u_perturb = u1
                        u_perturb(ivar) = u1(ivar) + perturbation_value
                        call perturbation_gradient(gradw1,gradw2, ivar, c1,c2, u_perturb, grad_perturb1,grad_perturb2)
                        
                        ! Compute the resultant changes in T and mu
                        T_perturb = w_perturb(5)*gamma/w_perturb(1)
                        mu_perturb = get_viscosity(T_perturb)
                        call perturb_temperature_gradient(gradT(:,c1),gradT(:,c2),c1,c2, T_perturb,gradT_perturb1,gradT_perturb2)
                        
                        call viscous_flux(w_perturb,w2, grad_perturb1,grad_perturb2, &
                                                                 T_perturb,Temp(c2), &
                                                gradT_perturb1(:),gradT_perturb2(:), &
                                                                  mu_perturb,mu(c2), &
                                                                   unit_face_normal, & !<- unit face normal
                                              cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                                              cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, & !<- Right cell centroid
                                                                     base_visc_flux  )
                        
                        dFnduL(:,ivar) = (num_flux - base_num_flux) / perturbation_value
                    else if(ii == 2) then

                    end if

                end do loop_var_perturbs_v
            end do loop_sides_v


        end do loop_faces




        ! Now we gotta do the boundary faces
        bound_loop : do ib = 1,nb
            bfaces_loop : do i = 1,bound(ib)%nbfaces
                c1 = bound(ib)%bcell(i)
                if (low_mach_correction) then
                    q1 = q(:,c1)
                    u1 = q2u(q1)
                    gradw1 = gradq(:,:,c1)
                else
                    u1 = u(1:5,c1)
                    gradw1 = gradw(1:3,1:5,c1)
                endif
                
                bface_centroid = bound(ib)%bface_center(:,i)
                unit_face_normal = bound(ib)%bface_nrml(:,i)
                face_mag       = bound(ib)%bface_nrml_mag(i)
                
                if (use_limiter) then
                    phi1 = phi(c1)
                    phi2 = one
                else 
                    phi1 = one
                    phi2 = one
                end if

                gradw2 = zero
                call get_right_state(bface_centroid(1), bface_centroid(2), bface_centroid(3), &
                                    u1, unit_face_normal, bc_type(ib), ub)
                wb = u2w(ub)
                if (low_mach_correction) then
                    qb = u2q(ub)
                    uR2b = min(qb(5) , max(qb(2)**2 + qb(3)**2 + qb(4)**2, epsilon*qb(5) ))
                    call interface_flux(          q1,      qb, & !<- Left/right states
                                        gradw1,        gradw2, & !<- Left/right gradients
                                             unit_face_normal, & !<- unit face normal
                        cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                                            bface_centroid(1), &
                                            bface_centroid(2), &
                                            bface_centroid(3), & !<- Right cell centroid
                                            bface_centroid(1), &
                                            bface_centroid(2), &
                                            bface_centroid(3), & !<- boundary ghost cell "center"
                                            phi1,        phi2, & !<- Limiter functions
                                    base_num_flux, wave_speed, & !<- Output)
                                       uR2L=uR2(c1), uR2R=uR2b )
                else
                    call interface_flux(          u1,      ub, & !<- Left/right states
                                        gradw1,        gradw2, & !<- Left/right gradients
                                             unit_face_normal, & !<- unit face normal
                        cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                                            bface_centroid(1), &
                                            bface_centroid(2), &
                                            bface_centroid(3), & !<- Right cell centroid
                                            bface_centroid(1), &
                                            bface_centroid(2), &
                                            bface_centroid(3), & !<- boundary ghost cell "center"
                                            phi1,        phi2, & !<- Limiter functions
                                    base_num_flux, wave_speed  ) !<- Output)
                endif
                loop_var_perturbs_b : do ivar = 1,5
                    ! u_perturb = u1
                    ! u_perturb(ivar) = u1(ivar) + perturbation_value
                    ! call perturbation_gradient_boundary(gradw1, ivar, c1, u_perturb, grad_perturb1)
                    ! call get_right_state(bface_centroid(1), bface_centroid(2), bface_centroid(3), &
                    !                 u_perturb, unit_face_normal, bc_type(ib), ub)
                    ! if (c1 == 1) then
                    !     write(*,*) c1
                    ! end if
                    if (low_mach_correction) then
                        call perturb_flux(   q1,      qb, & !<- Left/right states
                                            gradw1, gradw2, & !<- Left/right gradients
                                                unit_face_normal, & !<- unit face normal
                            cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                                                bface_centroid(1), &
                                                bface_centroid(2), &
                                                bface_centroid(3), & !<- Right cell centroid
                                                bface_centroid(1), &
                                                bface_centroid(2), &
                                                bface_centroid(3), & !<- boundary ghost cell "center"
                                                phi1,        phi2, & !<- Limiter functions
                                     ivar, perturbation_value,  1, &
                                             num_flux, wave_speed, & !<- Output)
                                            uR2L=uR2(c1), uR2R=uR2b )
                    else
                        call perturb_flux(   u1,      ub, & !<- Left/right states
                                          gradw1, gradw2, & !<- Left/right gradients
                                        unit_face_normal, & !<- unit face normal
                   cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                                       bface_centroid(1), &
                                       bface_centroid(2), &
                                       bface_centroid(3), & !<- Right cell centroid
                                       bface_centroid(1), &
                                       bface_centroid(2), &
                                       bface_centroid(3), & !<- boundary ghost cell "center"
                                       phi1,        phi2, & !<- Limiter functions
                            ivar, perturbation_value,  1, &
                                    num_flux, wave_speed  ) !<- Output)
                    endif
                    dFnduL(:,ivar) = (num_flux - base_num_flux) / perturbation_value

                end do loop_var_perturbs_b

                ! We only have a diagonal term to add
                jac(c1)%diag            = jac(c1)%diag            + dFnduL * face_mag

            end do bfaces_loop

        end do bound_loop

        
        !--------------------------------------------------------------------------------
        ! Add pseudo time term vol/dtau to the diagonals of the diagonal blocks
        ! to form: Jacobian = V/dtau + dRes1/dU, where V/dtau is a global diagonal
        ! matrix having vol(i)/dtau(i) for each node, and Res1 is the first-order
        ! residual.
        do i = 1,ncells
            if (low_mach_correction) then
                
                H = ((q(5,i))**2)/(gamma-one) + half * ( q(2,i)**2 + q(3,i)**2 + q(4,i)**2 )
                rho_p = gamma/q(5,i)
                rho_T = - (q(1,i)*gamma)/(q(5,i)**2)
                rho = q(1,i)*gamma/q(5,i)
                theta = (1/uR2(i)) - rho_T*(gamma - one)/(rho)
                
                preconditioner(1,:) = (/ theta,        zero,       zero,       zero,       rho_T                    /)
                preconditioner(2,:) = (/ theta*q(2,i), rho,        zero,       zero,       rho_T*q(2,i)             /)
                preconditioner(3,:) = (/ theta*q(3,i), zero,       rho,        zero,       rho_T*q(3,i)             /)
                preconditioner(4,:) = (/ theta*q(4,i), zero,       zero,       rho,        rho_T*q(4,i)             /)
                preconditioner(5,:) = (/ theta*H-one,  rho*q(2,i), rho*q(3,i), rho*q(4,i), rho_T*H + rho/(gamma-one)/)


                beta = rho_p + rho_T*(gamma-one)/rho
                alpha = half * (one - beta * uR2(i))
                absu = sqrt(q(2,i)**2 + q(3,i)**2 + q(4,i)**2 )
                lambda = absu*(one-alpha) + sqrt((alpha**2) * (absu**2) + uR2(i))
                dtau(i) = dtau(i) * two * wsn(i) / lambda

                jac(i)%diag(:,:) = jac(i)%diag(:,:) + preconditioner * cell(i)%vol/dtau(i)

            else
                do k = 1,5
                    jac(i)%diag(k,k) = jac(i)%diag(k,k) + cell(i)%vol/dtau(i)
                end do
            end if
            
        end do

        ! Invert the diagonal block matrices using Gauss elimination with pivoting
        !  Note: gewp_solve() actually solves a linear system, Ax=b, but here
        !        we use it only to obtain the inverse of A. So, b is a dummy.
        cell_inverse_loop : do i = 1,ncells
            !jac(i)%diag_inv = zero
            idestat = 0
            !                A                 dim  A^{-1}           error check
            call gewp_solve( jac(i)%diag(:,:), 5  , jac(i)%diag_inv, idestat    )
             !  Report errors
            if (idestat/=0) then
                write(*,*) " Error in inverting the diagonal block... Stop"
                write(*,*) "  Cell number = ", i
                do k = 1, 5
                    write(*,'(12(es8.1))') ( jac(i)%diag(k,j), j=1,5 )
                end do
                stop
            endif
        end do cell_inverse_loop

    end subroutine compute_numerical_jacobian

end module module_numerical_jacobian   