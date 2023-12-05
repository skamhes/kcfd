module module_jacobian
    implicit none
    
    public :: compute_jacobian ! compute the jacobian
    !public :: gewp_solve       ! Gauss elimination for inverting diagonal blocks
contains
    
    subroutine compute_jacobian
        
        use module_common_data          , only : p2, zero, nb, bc_type, half, one, two
        use module_flux_jac_interface   , only : interface_jac, interface_viscous_alpha_jacobian, & 
                                                 interface_flux_viscous_prim
        use module_ccfv_data_grid       , only : nfaces, face, ncells, cell, face_nrml, face_nrml_mag, jac, &
                                                 kth_nghbr_of_1, kth_nghbr_of_2, bound
        use module_ccfv_data_soln       , only : w, u, u2w, dtau, Temp, mu, gradw, uR2, gamma, wsn, q, u2q, res, q2u, mu, gradq
        use module_bc_states            , only : get_right_state, get_viscous_right_state
        use module_gewp                 , only : gewp_solve
        use module_input_parameter      , only : navier_stokes, low_mach_correction, visc_flux_method
        
        implicit none
        ! Local Vars
        integer                     :: c1, c2, i, k, ib, idestat, j, os
        real(p2), dimension(3)      :: unit_face_nrml, bface_centroid, ds, d_Cb, ejk
        real(p2), dimension(5)      :: u1, u2, ub, wb, qb, q1
        real(p2), dimension(3,5)    :: gradw1, gradw2, gradwb
        real(p2), dimension(5,5)    :: dFnduL, dFnduR
        real(p2)                    :: face_mag, mag_ds, mag_ejk
        real(p2)                    :: xc1,xc2,yc1,yc2,zc1,zc2

        real(p2), dimension(5,5)    :: preconditioner, pre_inv
        real(p2), dimension(5,5)    :: duLdqL, duRdqR
        real(p2)                    :: theta
        real(p2)                    :: rho_p, rho_T, rho
        real(p2)                    :: H, alpha, beta, lambda, absu
        
        ! Initialize jacobian terms
        do i = 1,ncells
            jac(i)%diag = zero
            jac(i)%off_diag = zero
            jac(i)%diag_inv = zero
        end do
        ! Loop Faces
        do i = 1,nfaces
            c1 = face(1,i)
            c2 = face(2,i)

            unit_face_nrml = face_nrml(1:3,i)
            face_mag       = face_nrml_mag(i)

            ! Compute the flux Jacobian for given w1 and w2
            if (low_mach_correction) then
                ! in this case dFndu_ is actually dFndq_.  I'm developing a rather long list of things that need to be fixed...
                call interface_jac( q(:,c1), q(:,c2), unit_face_nrml, dFnduL, dFnduR, uR2L = uR2(c1), uR2R = uR2(c2))
            else
                call interface_jac( w(:,c1), w(:,c2), unit_face_nrml, dFnduL, dFnduR)
            end if

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
            
            if (.not. navier_stokes) then
                cycle
            else
                if (visc_flux_method == 'alpha') then
                    u1      = u(:,c1)
                    u2      = u(:,c2)

                    gradw1  = gradw(:,:,c1)
                    gradw2  = gradw(:,:,c2)

                    ds      = (/cell(c2)%xc - cell(c1)%xc,cell(c2)%yc-cell(c1)%yc,cell(c2)%zc-cell(c1)%zc/)
                    mag_ds   = sqrt(ds(1)**2 + ds(2)**2 + ds(3)**2)
                    ds      = ds/mag_ds

                    call interface_viscous_alpha_jacobian(u1,u2,gradw1,gradw2,unit_face_nrml,ds,mag_ds,dFnduL,dFnduR)
                elseif (visc_flux_method == 'corrected') then
                    call interface_flux_viscous_prim(q(:,c1),q(:,c2),mu(c1),mu(c2),gradq(:,:,c1),gradq(:,:,c2),unit_face_nrml, &
                                                                                          cell(c1)%xc,cell(c1)%yc,cell(c1)%zc, & 
                                                                                          cell(c2)%xc,cell(c2)%yc,cell(c2)%zc, &
                                                                                          dFnduL,dFnduR)
                endif
                
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

            end if


        end do
        
        ! Now we gotta do the boundary faces
        bound_loop : do ib = 1,nb
            bfaces_loop : do i = 1,bound(ib)%nbfaces
                c1 = bound(ib)%bcell(i)
                
                bface_centroid = bound(ib)%bface_center(:,i)
                unit_face_nrml = bound(ib)%bface_nrml(:,i)
                face_mag       = bound(ib)%bface_nrml_mag(i)
                
                
                
                if (low_mach_correction) then
                    q1 = q(:,c1)
                    u1 = q2u(q1)
                    call get_right_state(bface_centroid(1), bface_centroid(2), bface_centroid(3), &
                                     u1, unit_face_nrml, bc_type(ib), ub)
                    qb = u2q(ub)
                    call interface_jac( q1, qb, unit_face_nrml, dFnduL, dFnduR, uR2L = uR2(c1), uR2R = uR2(c1))

                else
                    u1 = u(1:5,c1)
                    call get_right_state(bface_centroid(1), bface_centroid(2), bface_centroid(3), &
                                     u1, unit_face_nrml, bc_type(ib), ub)
                    wb = u2w(ub)
                    call interface_jac( w(:,c1), wb, unit_face_nrml, dFnduL, dFnduR)
                end if

                ! We only have a diagonal term to add
                jac(c1)%diag            = jac(c1)%diag            + dFnduL * face_mag

                if (.not. navier_stokes) then
                    cycle
                else
                    call get_viscous_right_state(bface_centroid(1),bface_centroid(2),bface_centroid(3), &
                    u1,gradw1,unit_face_nrml,bc_type(ib),ub,gradwb)

                    if (visc_flux_method == 'alpha') then
                        d_Cb = (/bface_centroid(1) - cell(c1)%xc,bface_centroid(2) - cell(c1)%yc,bface_centroid(3) - cell(c1)%zc/)
                        ejk = 2*d_Cb
                        mag_ejk = sqrt( ejk(1)**2 + ejk(2)**2 + ejk(3)**2 )
                        ejk = ejk/mag_ejk

    
                        call interface_viscous_alpha_jacobian(u1,ub,gradw1,gradwb,unit_face_nrml,ejk,mag_ejk,dFnduL,dFnduR)
                    elseif(visc_flux_method == 'corrected') then
                        call interface_flux_viscous_prim(q1,qb,mu(c1),mu(c1),gradq(:,:,c1),gradwb,unit_face_nrml, &
                        cell(c1)%xc,cell(c1)%yc,cell(c1)%zc, &
                        bface_centroid(1) - cell(c1)%xc, &
                        bface_centroid(2) - cell(c1)%yc, & 
                        bface_centroid(3) - cell(c1)%zc, &
                        dFnduL, dFnduR)
                    endif
                    ! if (low_mach_correction) then
    
                    !     ! We want dF/dQ = (dF/dU)(dU/dQ).  We don't need to redifine the left and right states.
                    !     dFnduL = matmul(dFnduL,duLdqL)

                    ! end if

                    jac(c1)%diag            = jac(c1)%diag            + dFnduL * face_mag

                end if
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
        
    end subroutine compute_jacobian


      
end module module_jacobian