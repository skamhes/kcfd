module module_ccfv_residual
    implicit none
    
    public :: compute_residual
contains
    subroutine compute_residual
        
        use module_common_data     , only : p2, zero, half, one

        use module_common_data     , only : x, y, bc_type, nb
       
        use module_ccfv_data_grid  , only : nfaces, face, ncells, cell, &
                    bound, face_nrml, face_nrml_mag, face_centroid
       
        use module_ccfv_data_soln  , only : res, u, gradw, wsn
       
        use module_flux            , only : interface_flux
        use module_bc_states       , only : get_right_state
        use module_ccfv_gradient   , only : compute_gradient
        use module_input_parameter , only : accuracy_order, use_limiter
        use module_ccfv_limiter    , only : compute_limiter, phi
        implicit none

        real(p2)                    :: xm, ym, zm
        integer                     :: i
        integer                     :: c1, c2

        real(p2), dimension(5)      :: u1, u2
        real(p2), dimension(3,5)    :: gradw1, gradw2
        real(p2), dimension(3)      :: unit_face_normal, bface_centroid

        real(p2), dimension(5)      :: num_flux
        integer                     :: j, ib, v1, v2, v3
        real(p2), dimension(5)      :: ub
        real(p2)                    :: wave_speed
        real(p2)                    :: phi1, phi2
        !real(p2), dimension(5,32)   :: local_res ! debugging
        !real(p2), dimension(3,5,32)   :: local_gradw ! debugging

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
        
        if (accuracy_order == 2) call compute_gradient
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
            !local_res = res
        end do loop_faces
        ! End of Residual computations for internal faces

        ! Now boundary faces
        boundary_loop : do ib = 1,nb
            bface_loop : do j = 1,bound(ib)%nbfaces
                bface_centroid = bound(ib)%bface_center(:,j)
                if (use_limiter) then
                    phi1 = phi(c1)
                    phi2 = phi(c2)
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
            end do bface_loop
        end do boundary_loop

    end subroutine compute_residual
end module module_ccfv_residual