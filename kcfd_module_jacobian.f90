module module_jacobian
    implicit none
    
    public :: compute_jacobian ! compute the jacobian
    !public :: gewp_solve       ! Gauss elimination for inverting diagonal blocks
contains
    
    subroutine compute_jacobian
        
        use module_common_data          , only : p2, zero, nb, bc_type
        use module_flux_jac_interface   , only : interface_jac, interface_jac_visc, interface_viscous_alpha_jacobian
        use module_ccfv_data_grid       , only : nfaces, face, ncells, cell, face_nrml, face_nrml_mag, jac, &
                                                 kth_nghbr_of_1, kth_nghbr_of_2, bound
        use module_ccfv_data_soln       , only : w, u, u2w, dtau, Temp, mu, gradw
        use module_bc_states            , only : get_right_state, get_viscous_right_state
        use module_gewp                 , only : gewp_solve
        use module_input_parameter      , only : navier_stokes
        
        implicit none
        ! Local Vars
        integer                     :: c1, c2, i, k, ib, idestat, j
        real(p2), dimension(3)      :: unit_face_nrml, bface_centroid, ds, d_Cb, ejk
        real(p2), dimension(5)      :: u1, u2, ub, wb
        real(p2), dimension(3,5)    :: gradw1, gradw2, gradwb
        real(p2), dimension(5,5)    :: dFnduL, dFnduR
        real(p2)                    :: face_mag, mag_ds, mag_ejk
        real(p2)                    :: xc1,xc2,yc1,yc2,zc1,zc2
        
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
            call interface_jac( w(:,c1), w(:,c2), unit_face_nrml, dFnduL, dFnduR)

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
                ! ! Calculate the face gradients
                ! xc1 = cell(c1)%xc
                ! yc1 = cell(c1)%yc
                ! zc1 = cell(c1)%zc
                ! xc2 = cell(c2)%xc
                ! yc2 = cell(c2)%yc
                ! zc2 = cell(c2)%zc
                ! delta_s = (/xc2-xc1, yc2-yc1, zc2-zc1/) ! vector pointing from center of cell 1 to cell 2

                ! call interface_jac_visc(w(:,c1),w(:,c2),Temp(c1),Temp(c2),mu(c1),mu(c2),gradw(:,2:4,c1),gradw(:,2:4,c2), &
                !                         unit_face_nrml,delta_s,dFnduL,dFnduR)

                ! ! Add to diagonal term of C1
                ! jac(c1)%diag            = jac(c1)%diag            + dFnduL * face_mag
                ! ! get neighbor index k for cell c1
                ! k = kth_nghbr_of_1(i)
                ! ! add to off diagonal neighbor k for cell c1
                ! jac(c1)%off_diag(:,:,k) = jac(c1)%off_diag(:,:,k) + dFnduR * face_mag

                ! ! Subtract terms from c2
                ! jac(c2)%diag            = jac(c2)%diag            - dFnduR * face_mag
                ! k = kth_nghbr_of_2(i)
                ! jac(c2)%off_diag(:,:,k) = jac(c2)%off_diag(:,:,k) - dFnduL * face_mag
                u1      = u(:,c1)
                u2      = u(:,c2)

                gradw1  = gradw(:,:,c1)
                gradw2  = gradw(:,:,c2)

                ds      = (/cell(c2)%xc - cell(c1)%xc,cell(c2)%yc-cell(c1)%yc,cell(c2)%zc-cell(c1)%zc/)
                mag_ds   = sqrt(ds(1)**2 + ds(2)**2 + ds(3)**2)
                ds      = ds/mag_ds

                call interface_viscous_alpha_jacobian(u1,u2,gradw1,gradw2,unit_face_nrml,ds,mag_ds,dFnduL,dFnduR)

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
                u1 = u(1:5,c1)
                bface_centroid = bound(ib)%bface_center(:,i)
                unit_face_nrml = bound(ib)%bface_nrml(:,i)
                face_mag       = bound(ib)%bface_nrml_mag(i)
                call get_right_state(bface_centroid(1), bface_centroid(2), bface_centroid(3), &
                                     u1, unit_face_nrml, bc_type(ib), ub)
                wb = u2w(ub)
                call interface_jac( w(:,c1), wb, unit_face_nrml, dFnduL, dFnduR)
                ! We only have a diagonal term to add
                jac(c1)%diag            = jac(c1)%diag            + dFnduL * face_mag

                if (.not. navier_stokes) then
                    cycle
                else

                    d_Cb = (/bface_centroid(1) - cell(c1)%xc,bface_centroid(2) - cell(c1)%yc,bface_centroid(3) - cell(c1)%zc/)
                    ejk = 2*d_Cb
                    mag_ejk = sqrt( ejk(1)**2 + ejk(2)**2 + ejk(3)**2 )
                    ejk = ejk/mag_ejk

                    call get_viscous_right_state(bface_centroid(1),bface_centroid(2),bface_centroid(3), &
                                                 u1,gradw1,unit_face_nrml,bc_type(ib),ub,gradwb)

                    call interface_viscous_alpha_jacobian(u1,ub,gradw1,gradwb,unit_face_nrml,ejk,mag_ejk,dFnduL,dFnduR)
                end if
            end do bfaces_loop
        end do bound_loop

        !--------------------------------------------------------------------------------
        ! Add pseudo time term vol/dtau to the diagonals of the diagonal blocks
        ! to form: Jacobian = V/dtau + dRes1/dU, where V/dtau is a global diagonal
        ! matrix having vol(i)/dtau(i) for each node, and Res1 is the first-order
        ! residual.
        do i = 1,ncells
            do k = 1,5
                jac(i)%diag(k,k) = jac(i)%diag(k,k) + cell(i)%vol/dtau(i)
            end do
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