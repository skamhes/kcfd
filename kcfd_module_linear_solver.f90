module module_linear_solver

    implicit none
    
    public :: linear_relaxation

contains

    subroutine linear_relaxation
        use module_common_data          , only : p2, lrelax_roc, lrelax_sweeps_actual, one, i_iteration, du, &
                                                 zero, my_eps, amg_direction, roc
        use module_input_parameter      , only : solver_type, lrelax_sweeps, lrelax_tolerance, lrelax_scheme
        use module_ccfv_data_soln       , only : res
        use module_ccfv_data_grid       , only : ncells, cell, jac

        implicit none
        
        real(p2), dimension(5)      :: b
        integer                     :: i, ii, k, os
        integer                     :: level = 1

        ! Linear system residual
        real(p2), dimension(5)      :: linear_res

        ! Residual norms(L1,L2,Linf)
        real(p2), dimension(5,3)    :: linear_res_norm, linear_res_norm_init, linear_res_norm_prev

        ! Rate of convergence
        real(p2) :: roc_previous

        ! Under-relaxation parameter
        real(p2) :: omega_lrelax

        ! Initialize some variables
        lrelax_sweeps_actual = 0 ! to count the actual number of relaxations
        omega_lrelax         = one
        amg_direction = 'up'
        do i = 1,ncells
            jac(i)%RHS = -res(:,i)
        end do

        call linear_sweeps(ncells,cell,jac,level,du)
        !Open the output file.
        !open(unit=1000, file='linear_solve.1000', status="unknown", iostat=os)
        ! Gonna skip writing the output for now...

        !if    ( trim(solver_type) == "implicit" ) then

        !    write(1000,*) "-------------------------------------------------------------------------------"
        !    write(1000,*) " Nonlinear Iteration = ", i_iteration
        !    write(1000,*) "-------------------------------------------------------------------------------"
           
        !elseif( trim(solver_type) == "jfnk" ) then
           
        !    write(1000,*) "-------------------------------------------------------------------------------"
        !    write(1000,*) " Nonlinear Iteration = ", i_iteration, " GCR projection = ", gcr_proj_actual
        !    write(1000,*) "-------------------------------------------------------------------------------"
           
        !endif
        !write(1000,*)
        !write(1000,*) " lrelax_scheme = ", trim(lrelax_scheme)
        !write(1000,*)

        

        lrelax_roc = roc
        !write(1000,*)
    
        return

        !close(1000)


    end subroutine linear_relaxation

    recursive subroutine linear_sweeps(ncells,cell,jac,level,correction)
        use module_common_data,     only : p2, zero, one, lrelax_sweeps_actual, my_eps, half, amg_direction, roc, two
        use module_ccfv_data_grid,  only : cc_data_type, jacobian_type
        use module_input_parameter, only : solver_type, lrelax_sweeps, lrelax_tolerance, lrelax_scheme, &
                                            use_amg, max_amg_levels

        implicit none
        integer, intent(in)                                 :: ncells
        type(cc_data_type), dimension(ncells), intent(in)   :: cell
        type(jacobian_type), dimension(ncells), intent(in)  :: jac
        integer, intent(in)                                 :: level
        real(p2), dimension(5,ncells), intent(out)          :: correction
        
        real(p2), dimension(5,ncells)                       :: additive_correction
        ! Residual norms(L1,L2,Linf)
        real(p2), dimension(5,3)    :: linear_res_norm, linear_res_norm_init, linear_res_norm_prev

        real(p2), dimension(5)      :: b
        integer                     :: i, ii, k, os
        
        ! Linear system residual
        real(p2), dimension(5)      :: linear_res

        ! Under-relaxation parameter
        real(p2) :: omega_lrelax

        ! Rate of convergence
        real(p2) :: roc_previous

        ! Initialize some variables
        omega_lrelax         = one

        !--------------------------------------------------------------------------------
        ! 1. Initialize the correction
    
        correction = zero

        !--------------------------------------------------------------------------------
        ! 2. Linear Relaxation (Sweep)
        relax : do ii = 1,lrelax_sweeps
            linear_res_norm(:,1) = zero
            !---------------------------------------------------------
            ! Sequential Gauss-Seidel Relaxation(sweep)
            if (trim(lrelax_scheme) == 'gs') then
                gs_loop : do i = 1,ncells
                    ! Form the right hand side of GS: [ sum( off_diagonal_block*du ) - residual ]
                    b = jac(i)%RHS
                    gs_nghbr_loop : do k = 1,cell(i)%nnghbrs
                        ! Add RHS from off diagonal terms and du (du = zero to start and will be updated as we go)
                        b = b - matmul(jac(i)%off_diag(:,:,k), correction(:,cell(i)%nghbr(k)))
                    end do gs_nghbr_loop 
                    ! Update du by the GS relaxation:
                    !
                    ! e.g., for 3 nghbrs, perform the relaxation in the form:
                    !
                    !                     diagonal block        sum of off-diagonal block contributions
                    !       dUj = omega*{ [V/dtj+dR/dUj]^{-1}*(-[dRj/dU1]*dU1 -[dRj/dU2]*dU2 -[dRj/dU3]*dU3 -Res_j) - dUj }
                    linear_res = matmul(jac(i)%diag_inv(:,:), b) - correction(:,i)
                    correction(:,i) = correction(:,i) + omega_lrelax * linear_res
                    linear_res_norm(:,1) = linear_res_norm(:,1) + abs(linear_res)
                end do gs_loop
            else
                write(*,*) " Sorry, only 'gs' is available at the moment..."
                write(*,*) " Set lrelax_scheme = 'gs', and try again. Stop."
                stop
            end if
            !---------------------------------------------------------
            linear_res_norm(:,1) = linear_res_norm(:,1) / real(ncells, p2)

            !--------------------------------------------------------------------------------
            ! 3. Check the linear residual.

            !  After the first relaxation
            if (ii == 1) then
                !-----------------------------------------------------------------
                ! Print the initial linear residual norm in the file 'fort.1000'.

                linear_res_norm_init = linear_res_norm
                !write(1000,'(a,i10,a,es12.5)') "  after relax ", ii, &
                !            " max(L1 norm) = ", maxval(linear_res_norm(:,1))

                if ( maxval(linear_res_norm(:,1)) < my_eps ) then
                    !write(1000,*) " Machine zero res reached. Exit GS relaxation. Total sweeps = ", ii
                    !write(1000,*) " tolerance_linear = ", lrelax_tolerance
                    lrelax_sweeps_actual = ii
                    exit relax
                endif
            ! After the second Relaxation
            else 
                roc = maxval(linear_res_norm(1:5,1)/linear_res_norm_init(1:5,1))

                !-----------------------------------------------------------------
                ! Print the current linear residual norm in the file 'fort.1000'.
                if (roc < one) then
                    !write(1000,'(a,i10,a,es12.5,2x,f6.3,2x,f8.3)')   "  after relax ", ii,          &
                    !            " max(L1 norm) = ", maxval(linear_res_norm(1:5,1)), omega_lrelax, roc
                else

                    !write(1000,'(a,i10,a,es12.5,2x,f6.3,2x,f8.3,a)') "  after relax ", ii,          &
                    !        " max(L1 norm) = ", maxval(linear_res_norm(1:5,1)), omega_lrelax, roc," <- diverge"
                endif
                !-----------------------------------------------------------------
                ! Tolerance met: Exit

                if (roc < lrelax_tolerance) then
                ! if (roc < lrelax_tolerance*two*level) then
                ! if (roc < min(lrelax_tolerance*(10**(level-1)),two)) then
                    !write(1000,*)
                    !write(1000,*) " Tolerance met. Exit GS relaxation. Total sweeps = ", ii
                    !write(1000,*) " tolerance_linear = ", lrelax_tolerance
                    lrelax_sweeps_actual = ii
                    ! write(*,*) ii
                    exit relax
                !-----------------------------------------------------------------
                ! If tolerance is NOT met.
                else

                    !---------------------------------------------
                    ! Stop if the linear residual is too small.
            
                    if ( maxval(linear_res_norm(1:5,1)) < my_eps ) then
                    !    write(1000,*)
                    !    write(1000,*) " Residuals too small. Exit GS relaxation. Total sweeps = ", ii
                    !    write(1000,*) " maxval(linear_res_norm(1:5,1)) = ", maxval(linear_res_norm(1:5,1))
                        lrelax_sweeps_actual = ii
                        exit relax
                    endif
            
                    !-------------------------------------------------
                    ! Stop if the maximum number of sweeps is reached.
            
                    if (ii == lrelax_sweeps) then
                    !    write(1000,*)
                    !    write(1000,*) " Tolerance not met... sweeps = ", lrelax_sweeps
                        lrelax_sweeps_actual = lrelax_sweeps
                    endif
                end if

                ! End of Tolerance met or NOT met.
                !-----------------------------------------------------------------

                !-----------------------------------------------------------------
                ! Do something if we continue.
                ! Reduce/increase the relaxation factor if it's going worse/better.

                roc_previous = maxval(linear_res_norm_prev(1:5,1)/linear_res_norm_init(1:5,1))

                !--------------------------------------------------------------
                ! If the linear residual goes beyond the initial one,
                if (roc > one) then
                    ! REDUCE if the linear residual increases from the previous relaxation.
                    if (roc > roc_previous) omega_lrelax = max(0.95_p2*omega_lrelax, 0.05_p2)
                    !--------------------------------------------------------------
                    ! If the linear residual is smaller than the initial one,
                else
                    ! INCREASE if doing fine.
                    omega_lrelax = min(1.05_p2*omega_lrelax, one)
                    if (roc < roc_previous) omega_lrelax = min(1.25_p2*omega_lrelax, one)    
                endif
            
                ! End of Do something if we continue.
                !-----------------------------------------------------------------
            
            endif
            linear_res_norm_prev = linear_res_norm
            
            if (use_amg .and. (level <= max_amg_levels) .and.  & 
                (ii > 10) .and. (roc/roc_previous > half) .and. (amg_direction == 'up')) then ! will probably play with the convergence condition here some
                call algebraic_multigrid(ncells,cell,correction, jac, level, additive_correction)
                do i = 1,ncells
                    correction(:,i) = correction(:,i) + additive_correction(:,i)
                end do
            end if
        end do relax
        ! End of 3. Linear Relaxation (Sweep)
        !--------------------------------------------------------------------------------

    end subroutine linear_sweeps
    
    recursive subroutine algebraic_multigrid(ncells,cell,phi,jac,level,correction_del)
        use module_common_data   , only : p2, zero, amg_direction
        use module_ccfv_data_grid, only : cc_data_type, jacobian_type, test
        use module_jacobian      , only : gewp_solve
        
        implicit none
        ! Input
        integer                                              :: ncells ! number of cells (level 1) or groups (level 2+)
        type(cc_data_type), dimension(ncells), intent(in)    :: cell   ! cell data (level 1) or group data stored as a cell (l2)
        !real(p2), dimension(ncells), intent(in)              :: dtau   ! dtau ! pretty sure I don't need this
        real(p2), dimension(5,ncells), intent(in)            :: phi    ! dependent variables (du for level 1, del for level 2+)
        type(jacobian_type), dimension(ncells), intent(in)   :: jac    ! jacobian matrix (level 1) restricted operator (level 2+)
        integer, intent(in)                                  :: level  ! Multigrid level (not sure if we'll need it yet)
        ! Output
        real(p2), dimension(5,ncells), intent(out)           :: correction_del ! prolongated correction factor (del) to be 
                                                                           ! passed to finer level

        real(p2), dimension(:,:), allocatable :: g_correction ! restricted correction factor
        real(p2)                            :: strongest_A11, group_size
        integer, dimension(ncells)          :: assign_group
        integer                             :: i, k, ngroup, ck, gi, ci, kk, idestat, strongest_k, last_k
        type(cc_data_type),  dimension(:), allocatable :: group_cell
        type(jacobian_type), dimension(:), allocatable :: group
        ! real(p2), dimension(5,5) :: current_diag, current_diag_inv
        ! real(p2),dimension(5,5,4) :: current_off_diag
        ! real(p2), dimension(5) :: current_rhs
        integer, dimension(1000)                        :: temp_nghbr, temp_gnghbr
        ! Initialize var
        assign_group = zero
        correction_del = zero
        ngroup = 0
        assign_cells : do i = 1,ncells
            if ((assign_group(i) /= 0)) cycle
            ! strongest_A11 = zero
            ! strongest_k = 0
            ngroup = ngroup + 1
            group_size = 1
            last_k = 0
            assign_group(i) = ngroup
            do k = 1,cell(i)%nnghbrs
                ! if ((abs(jac(i)%off_diag(1,1,k)) > abs(strongest_A11)) .and. (assign_group(cell(i)%nghbr(k)) == 0)) then
                !     strongest_A11 =  jac(i)%off_diag(1,1,k)
                !     strongest_k = k
                ! end if
                if (assign_group(cell(i)%nghbr(k)) == 0) then
                    assign_group(cell(i)%nghbr(k)) = ngroup
                    group_size = group_size + 1
                    last_k = k
                end if
                if (group_size >= 4) cycle assign_cells
            end do
            if ( last_k == 0 .and. cell(i)%nnghbrs /= 0) then ! all neighbors are in a group
                assign_group(i) = assign_group(cell(i)%nghbr(1)) ! just add it to the neighbor of group 1 (theoretically this 
                ! could result in groups as large as 6 but this is unlikely and the cut off of 4 is somewhat arbitrary)
                ngroup = ngroup - 1 ! remove the group
                cycle assign_cells
            end if
            ! we only make it here if group_size < 4
            ! so we're just gonna keep adding cells till we get to 4
            do kk = 1,cell(cell(i)%nghbr(last_k))%nnghbrs
                if (assign_group(cell(cell(i)%nghbr(last_k))%nghbr(kk)) == 0) then
                    assign_group(cell(cell(i)%nghbr(last_k))%nghbr(kk)) = ngroup
                    group_size = group_size + 1
                end if
                if (group_size >= 4) cycle assign_cells
            end do
            ! if (strongest_k == 0) then
            !     ngroup = ngroup + 1
            !     assign_group(i) = ngroup ! group of 1
            ! else 
            !     ngroup = ngroup + 1
            !     assign_group(i) = ngroup
            !     assign_group(cell(i)%nghbr(strongest_k)) = ngroup
            ! end if
        end do assign_cells
        ! Allocate group arrays
        allocate(group(ngroup))
        allocate(group_cell(ngroup))
        ! We're going to take advantage of the unused ntvx and vtx vars to store the # of cells and cells in each group
        ! This is fine as it's only being used internally (don't have to worry about consistency when passing between levels)
        do i = 1,ngroup
            group_cell(i)%nvtx = 0
        end do
        do i = 1,ncells
            gi = assign_group(i)
            group_cell(gi)%nvtx = group_cell(gi)%nvtx + 1
        end do
        do gi = 1,ngroup
            allocate(group_cell(gi)%vtx(group_cell(gi)%nvtx))
            group_cell(gi)%vtx = 0
        end do
        do i = 1,ncells
            gi = assign_group(i)
            k = 1
            add_vtx : do
                if (group_cell(gi)%vtx(k) == 0) then
                    group_cell(gi)%vtx(k) = i    
                    exit add_vtx
                end if
                k = k + 1
            end do add_vtx
            ! if (group_cell(gi)%vtx(1) == 0) then
            !     group_cell(gi)%vtx(1) = i
            ! else 
            !     ! if (gi == 182 .and. level == 3) then
            !     !     write(*,*) gi, level
            !     ! end if
            !     group_cell(gi)%vtx(2) = i ! max group size is 2 so this is fine (for now) not anymore
            ! end if
        end do
        do gi = 1,ngroup
            group_cell(gi)%nnghbrs = 0
            group_cell(gi)%cnnghbrs = 0
            temp_nghbr = 0
            temp_gnghbr = 0
            do k = 1,group_cell(gi)%nvtx
                ck = group_cell(gi)%vtx(k)
                do kk = 1,cell(ck)%nnghbrs
                    if ((assign_group(cell(ck)%nghbr(kk)) /= gi) ) then
                        group_cell(gi)%nnghbrs = group_cell(gi)%nnghbrs + 1
                        temp_nghbr(group_cell(gi)%nnghbrs) = cell(ck)%nghbr(kk)
                        if (any(temp_nghbr /= cell(ck)%nghbr(kk))) then
                            group_cell(gi)%cnnghbrs = group_cell(gi)%cnnghbrs + 1
                            temp_gnghbr(group_cell(gi)%cnnghbrs) = assign_group(cell(ck)%nghbr(kk))
                        end if
                    end if
                end do
            end do
            allocate(group_cell(gi)%nghbr( group_cell(gi)%nnghbrs ))
            allocate(group_cell(gi)%cnghbr(group_cell(gi)%cnnghbrs))
            group_cell(gi)%nghbr  = temp_gnghbr( 1:group_cell(gi)%nnghbrs )
            group_cell(gi)%cnghbr = temp_nghbr(1:group_cell(gi)%cnnghbrs)
        
            ! Initialize group matrix
            group(gi)%diag = zero
            allocate(group(gi)%off_diag(5,5,group_cell(gi)%nnghbrs))
            group(gi)%off_diag = zero
            group(gi)%rhs = zero
        ! end do

        ! Add terms to restricted matrix (group)
        ! do i = 1,ncells
        !     gi = assign_group(i)
        !     ! group(gi)%diag = group(gi)%diag + jac(i)%diag
        !     ! group(gi)%RHS  = group(gi)%RHS  - matmul(jac(i)%diag,phi(:,i)) + jac(i)%RHS
        !     do k = 1,cell(i)%nnghbrs
        !         ! if (assign_group(cell(i)%nghbr(k)) == gi) then
        !         !     group(gi)%diag = group(gi)%diag + jac(i)%off_diag(:,:,k)
        !         ! end if
        !         ! group(gi)%RHS = group(gi)%RHS - matmul(jac(i)%off_diag(:,:,k),phi(:,i))
        !     end do
        ! end do
        ! Loop through the groups to add one more term
        ! do gi = 1,ngroup
            ! current_off_diag = zero
            i = 0
            do k = 1,group_cell(gi)%nvtx
                ck = group_cell(gi)%vtx(k)
                group(gi)%diag = group(gi)%diag + jac(ck)%diag
                group(gi)%RHS  = group(gi)%RHS  - matmul(jac(ck)%diag,phi(:,ck)) + jac(ck)%RHS
                do kk = 1,cell(ck)%nnghbrs
                    if (.not.(assign_group(cell(ck)%nghbr(kk)) == gi)) then
                        i = i + 1
                        ! write(*,*) gi, i
                        group(gi)%off_diag(:,:,i) = group(gi)%off_diag(:,:,i) + jac(ck)%off_diag(:,:,kk)
                        ! current_off_diag(:,:,i) = group(gi)%off_diag(:,:,i) 
                    end if
                    if (assign_group(cell(ck)%nghbr(kk)) == gi) then
                        group(gi)%diag = group(gi)%diag + jac(ck)%off_diag(:,:,kk)
                    end if
                    group(gi)%RHS = group(gi)%RHS - matmul(jac(ck)%off_diag(:,:,kk),phi(:,ck))
                end do
            end do
            call gewp_solve( group(gi)%diag(:,:), 5  , group(gi)%diag_inv, idestat    )
            ! current_diag_inv = group(gi)%diag_inv
        end do

        ! Call the linear solver to relax the newly coarsened system.
        allocate(g_correction(5,ngroup))
        call linear_sweeps(ngroup, group_cell, group, level + 1, g_correction)
        
        ! Add correction to phi
        do gi = 1,ngroup
            do k = 1,group_cell(gi)%nvtx
                ck = group_cell(gi)%vtx(k)
                correction_del(:,ck) = correction_del(:,ck) + g_correction(:,gi)
            end do
            deallocate(group_cell(gi)%vtx)
            deallocate(group_cell(gi)%nghbr)
            deallocate(group_cell(gi)%cnghbr)
            deallocate(group(gi)%off_diag)
        end do

        deallocate(group)
        deallocate(group_cell)
        deallocate(g_correction)
        ! Deallocate
        !call destroy_cells(ngroup, group)
        amg_direction = 'dn' ! we've finished a multigrid level that means we're going back down the levels
    end subroutine algebraic_multigrid
end module module_linear_solver
