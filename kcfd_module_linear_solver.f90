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
        use module_sparse_matrix        , only : build_A_BCSM

        implicit none
        
        real(p2), dimension(:,:,:), allocatable :: V   ! Values (5x5 block matrix) plus corresponding index
        integer , dimension(:),     allocatable :: C   ! Column index of each value
        integer , dimension(ncells+1)           :: R   ! Start index of each new row
        integer                                 :: nnz
        real(p2), dimension(5,5,ncells)             :: Dinv

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

        call build_A_BCSM(ncells,cell,jac,V,C,R,nnz=nnz)

        call build_Dinv_array(ncells,jac,Dinv)
        
        call linear_sweeps(ncells,V,C,R,nnz,res,Dinv,level,du)
        ! call linear_sweeps(ncells,cell,V,C,R,nnz,res,Dinv,level,du)
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
        
        ! We have a bunch of stuff to deallocate now...
        if (allocated(V)) deallocate(V)
        if (allocated(C)) deallocate(C)

        return

        !close(1000)


    end subroutine linear_relaxation

    recursive subroutine linear_sweeps(ncells,V,C,R,nnz,res,Dinv,level,correction)
        ! This subroutine relaxes the system using a Gauss-Seidel method with Algebraic Multigrid.

        use module_common_data,     only : p2, zero, one, lrelax_sweeps_actual, my_eps, half, amg_direction, roc, two, fourth
        use module_ccfv_data_grid,  only : cc_data_type, jacobian_type
        use module_input_parameter, only : solver_type, lrelax_sweeps, lrelax_tolerance, lrelax_scheme, &
                                            use_amg, max_amg_levels
        use module_sparse_matrix,   only : sparseblock_times_vectorblock, build_A_BCSM

        implicit none
        integer, intent(in)                                 :: ncells
        ! type(cc_data_type), dimension(ncells), intent(in)   :: cell
        ! type(jacobian_type), dimension(ncells), intent(in)  :: jac
        real(p2), dimension(:,:,:),allocatable, intent(in)  :: V    ! Values of A
        integer , dimension(:), allocatable, intent(in)     :: C    ! Column index of A
        integer , dimension(ncells+1), intent(in)           :: R    ! Start index of A
        integer , intent(in)                                :: nnz  ! # non-zero terms in A
        real(p2), dimension(5,ncells), intent(in)           :: res  ! RHS (= -b)
        real(p2), dimension(5,5,ncells), intent(in)         :: Dinv ! Inverse of A(i,i)
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
                gs_loop : do i = 1,ncells ! loop through rows
                    ! Form the right hand side of GS: [ sum( off_diagonal_block*du ) - residual ]
                    b = -res(:,i)
                    gs_row_loop : do k = R(i),(R(i+1)-1)
                        ! Add RHS from off diagonal terms and du (du = zero to start and will be updated as we go)
                        if (C(k) /= i) then
                            b = b - matmul(V(:,:,k),correction(:,C(k)))
                        end if
                    end do gs_row_loop
                    ! ! Update du by the GS relaxation:
                    !
                    ! e.g., for 3 nghbrs, perform the relaxation in the form:
                    !
                    !                     diagonal block        sum of off-diagonal block contributions
                    !       dUj = omega*{ [V/dtj+dR/dUj]^{-1}*(-[dRj/dU1]*dU1 -[dRj/dU2]*dU2 -[dRj/dU3]*dU3 -Res_j) - dUj }
                    linear_res = matmul(Dinv(:,:,i), b) - correction(:,i)
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

                ! if (roc < lrelax_tolerance) then
                if (roc < lrelax_tolerance/(2.25_p2*level-1.25_p2)) then
                                            ! level = 1 => lrelax, level = 5 => lrelax/10
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
                (ii > 5) .and. (roc/roc_previous > half) .and. (amg_direction == 'up')) then 
                call algebraic_multigrid(ncells,correction, V, C, R, nnz, res, level, additive_correction)
                do i = 1,ncells
                    correction(:,i) = correction(:,i) + additive_correction(:,i)
                end do
            end if
        end do relax
        ! End of 3. Linear Relaxation (Sweep)
        !--------------------------------------------------------------------------------
        ! write(*,*) ii
    end subroutine linear_sweeps
    
    recursive subroutine algebraic_multigrid(ncells,phi,V,C,R,nnz,res,level,correction_del)
    
        ! Implementation of algebraic multi grid using additive correction.
        use module_common_data   , only : p2, zero, amg_direction
        use module_ccfv_data_grid, only : cc_data_type, jacobian_type, test, jac
        use module_jacobian      , only : gewp_solve
        use module_sparse_matrix , only : R_A_P
        !use stdlib_sorting, only : sort
        
        implicit none
        ! Input
        integer                                              :: ncells ! number of cells (level 1) or groups (level 2+)
        real(p2), dimension(5,ncells), intent(in)            :: phi    ! dependent variables (du for level 1, del for level 2+)
        real(p2), dimension(:,:,:), intent(in)               :: V   ! Values (5x5 block matrix) plus corresponding index
        integer, dimension(:), intent(in)                    :: C   ! Column index of each value
        integer,dimension(ncells+1), intent(in)              :: R   ! Start index of each new row
        integer, intent(in)                                  :: nnz
        real(p2), dimension(5,ncells), intent(in)            :: res    ! RHS of the equation
        integer, intent(in)                                  :: level  ! Multigrid level (not sure if we'll need it yet)
        
        ! Output
        real(p2), dimension(5,ncells), intent(out)           :: correction_del ! prolongated correction factor (del) to be 
                                                                           ! passed to finer level

        ! Local Vars
        real(p2), dimension(:,:), allocatable :: g_correction ! restricted correction factor
        real(p2)                            :: strongest_A11, group_size, current_n
        real(p2), dimension(3)              :: A11_nghbrs
        integer,  dimension(3)              :: strongest_3nghbr, sort_index
        integer, dimension(ncells)          :: assign_group
        integer                             :: i, j, k, ck, gi, kk, os, ii, nz_count
        integer                             :: ngroup
        integer                             :: idestat

        real(p2), dimension(5,ncells)       :: defect ! defect used in RHS for AMG
        real(p2), dimension(:,:), allocatable :: defect_res ! R(A*phi + b)

        ! Restricted Vars
        integer, dimension(ncells)              :: RestrictC, ProlongC
        integer, dimension(ncells + 1)          :: RestrictR, ProlongR ! Restrict array is longer than needed but this proves fine.
        real(p2), dimension(:,:,:), allocatable :: RAP_V
        integer, dimension(:), allocatable      :: RAP_C
        integer, dimension(:), allocatable      :: RAP_R
        real(p2), dimension(:,:,:), allocatable :: RAP_Dinv

        ! Initialize var
        assign_group = 0
        correction_del = zero
        ngroup = 0
        RestrictR(1) = 1
        ProlongR = (/ (i, i = 1,(ncells+1)) /)
        assign_cells : do i = 1,ncells
            if ((assign_group(i) /= 0)) cycle
            ! Initialize
            ngroup = ngroup + 1
            group_size = 1
            A11_nghbrs = zero
            strongest_3nghbr = 0
            assign_group(i) = ngroup
            do j = R(i),(R(i+1)-1)
                if (C(j) == i .or. assign_group(C(j)) /= 0) cycle ! skip diagonal and already assigned
                current_n = ABS(V(1,1,j))
                if (current_n > A11_nghbrs(1)) then
                    A11_nghbrs(1) = current_n
                    strongest_3nghbr(1) = C(j)
                    call insertion_sort_real(A11_nghbrs,sort_index)
                    A11_nghbrs(:) = A11_nghbrs(sort_index)
                    strongest_3nghbr(:) = strongest_3nghbr(sort_index)
                end if
            end do
            if ( all(strongest_3nghbr(:) == 0)) then ! all neighbors are in a group
                ! Were just making a group of 1.  May look to add to a neighboring group in the future.  Adds some complexity to the
                ! CSRM foramtion of the restriction matrix.
                ProlongC(i) = ngroup
                RestrictR(ngroup + 1)        = RestrictR(ngroup) + 1 
                RestrictC(RestrictR(ngroup)) = i
                cycle assign_cells
            end if
            ! We only make it here if group_size < 4
            ! So we're just going to keep adding cells till we get to 4
            do j = 1,3
                if ( strongest_3nghbr(j) /= 0 ) then
                    assign_group(strongest_3nghbr(j)) = ngroup
                end if
            end do
            if ( all(strongest_3nghbr(:) /= 0) ) then
                ProlongC(i) = ngroup
                do j = 1,3
                    ProlongC(strongest_3nghbr(j)) = ngroup
                end do
                RestrictR(ngroup + 1) = RestrictR(ngroup) + 4
                RestrictC(RestrictR(ngroup):(RestrictR(ngroup+1)-1)) = insertion_sort_int((/i,strongest_3nghbr/),4)
                cycle ! if all 3 neighbors have been assigned we don't need to look for more
                ! only get here if they haven't been
            end if
            nghbr_group_add : do j = R(strongest_3nghbr(3)) , (R(strongest_3nghbr(3) + 1) - 1) 
                ! move to row of strongest neighbor
                if ( assign_group(C(j)) /= 0 .or. j == C(R(i)) ) cycle ! allready assigned or diagonal
                current_n = ABS(V(1,1,j))
                if (current_n > A11_nghbrs(1)) then
                    A11_nghbrs(1) = current_n
                    strongest_3nghbr(1) = C(j)
                    call insertion_sort_real(A11_nghbrs,sort_index)
                    A11_nghbrs(:) = A11_nghbrs(sort_index)
                    strongest_3nghbr(:) = strongest_3nghbr(sort_index)
                end if
                if ( all(strongest_3nghbr(:) /= 0)) exit nghbr_group_add ! only adding neighbors, not replacing.  Once we fill up
                ! were done.
            end do nghbr_group_add 
            ! add any missing terms
            nz_count = 0
            do j = 1,3
                if ( strongest_3nghbr(j) /= 0 ) then
                    assign_group(strongest_3nghbr(j)) = ngroup
                    ProlongC(strongest_3nghbr(j)) = ngroup
                    nz_count = nz_count + 1
                end if
            end do
            ProlongC(i) = ngroup
            RestrictR(ngroup + 1) = RestrictR(ngroup) + nz_count + 1
            if (nz_count > 0) then
                RestrictC(RestrictR(ngroup):(RestrictR(ngroup+1)-1)) = &
                    insertion_sort_int((/i,strongest_3nghbr((4-nz_count):)/),nz_count + 1)
            else
                RestrictC(RestrictR(ngroup)) = i ! I don't think it's possible to get this far with nz_count = 0 but I'm leaving
                ! this here just in case...
            end if
        end do assign_cells

        ! Create coarse level operetor A^H = RAP
        allocate(RAP_R(ngroup + 1))
        call R_A_P(ncells,ngroup,nnz,RestrictC,RestrictR,ProlongC,ProlongR,V,C,R,RAP_V,RAP_C,RAP_R,os)
        
        ! Create inverse block matrix of RAP diagonal terms
        allocate(RAP_Dinv(5,5,ngroup))
        do gi = 1,ngroup
            do j = RAP_R(gi),(RAP_R(gi+1)-1)
                if (RAP_C(j) == gi) then ! diag
                    call gewp_solve(RAP_V(:,:,j),5, RAP_Dinv(:,:,gi),idestat)
                end if
            end do
        end do
        
        ! Calculate and restrict the defect d = A*phi + b
        call compute_defect(ncells,V,C,R,phi,res,defect)
        allocate(defect_res(5,ngroup))
        defect_res = zero
        do i = 1,ngroup
            do j = RestrictR(i),(RestrictR(i+1)-1)
                defect_res(:,i) = defect_res(:,i) + defect(:,RestrictC(j))
            end do
        end do

        ! Call the linear solver to relax the newly coarsened system and apply the correction.
        allocate(g_correction(5,ngroup))
        call linear_sweeps(ngroup, RAP_V,RAP_C,RAP_R,nnz,defect_res,RAP_Dinv, level + 1, g_correction)
        do i = 1,ncells
            ! Since ProlongC has 1 value per row we can skip the inner j loop.
            correction_del(:,i) = g_correction(:,ProlongC(i))
        end do

        ! Make sure all of the allocated arrays are deallocated
        if (allocated(RAP_V))        deallocate(     RAP_V)
        if (allocated(RAP_C))        deallocate(     RAP_C)
        if (allocated(RAP_R))        deallocate(     RAP_R)
        if (allocated(RAP_Dinv))     deallocate(  RAP_Dinv)
        if (allocated(defect_res))   deallocate(defect_res)
        if (allocated(g_correction)) deallocate(g_correction)

        amg_direction = 'dn' ! we've finished a multigrid level that means we're going back down the levels
    end subroutine algebraic_multigrid

    subroutine create_restriction_matrix(ncells,ngroups,assign_group,Restriction)
        ! Builds a rank 2 restriction matrix by recognizing that the row of R defines the group and the column defines the cell
        ! Currently not used by kept around for debugging purposes in case it is needed.
        use  module_common_data , only : p2

        implicit none

        integer, intent(in) :: ncells
        integer, intent(in) :: ngroups
        integer, dimension(ncells) :: assign_group

        integer, dimension(ngroups,ncells), INTENT(OUT) :: Restriction

        integer :: i
        restriction = 0
        do i = 1,ncells
            restriction(assign_group(i),i) = 1
        end do

    end subroutine create_restriction_matrix

    subroutine compute_defect(ncells,V,C,R,phi,b,defect)
        ! Computes the defect d = A*phi + b
        use module_common_data , only : p2
        use module_sparse_matrix , only : sparseblock_times_vectorblock

        implicit none

        integer, intent(in) :: ncells
        real(p2), dimension(:,:,:), intent(in) :: V
        integer, dimension(:), intent(in) :: C
        integer, dimension(ncells + 1) :: R
        real(p2), dimension(5,ncells), intent(in) :: phi
        real(p2), dimension(5,ncells), intent(in) :: b

        real(p2), dimension(5,ncells), intent(out) :: defect

        real(p2), dimension(5,ncells) :: product

        call sparseblock_times_vectorblock(ncells,V,C,R,phi,product)

        defect = product + b

    end subroutine compute_defect

    subroutine restrict_defect(ncells,ngroups,defect,restrict,product)
        ! Compute the restricted defect R*d
        ! Note: this function is not currently used.  Kept for debugging purposes (same as create_restriction_matrix subroutine)
        use module_common_data , only : p2

        implicit none
        integer, intent(in) :: ncells
        integer, intent(in) :: ngroups
        real(p2), dimension(5,ncells), intent(in) :: defect
        real(p2), dimension(ngroups,ncells), intent(in) :: restrict

        real(p2), dimension(5,ngroups), intent(out) :: product

        integer :: i

        do i = 1,5
            product(i,:) = matmul(restrict,defect(i,:)) ! whether this works will depend on if matmul treats defect(i,:) as a rank
            ! one array.  We will see in debugging...
        end do
    end subroutine restrict_defect

    subroutine build_Dinv_array(ncells,jac,D_inv)
        ! This subroutine stores just the diagonal blocks of the Dinv matrix since the rest are empty.  As a result, C=R=index
        ! so the other two indices do not need to be stored.
        use module_common_data , only : p2, zero
        use module_ccfv_data_grid , only : jacobian_type

        implicit none
        integer, intent(in) :: ncells
        type(jacobian_type), dimension(ncells), intent(in) :: jac
        
        real(p2), dimension(5,5,ncells), INTENT(OUT) :: D_inv
        

        integer :: i
        
        do i = 1,ncells
            D_inv(:,:,i) = jac(i)%diag_inv(:,:)
        end do
    end subroutine build_Dinv_array

    subroutine LU_extraction(ncells,nnz,V,C,R,LU_V,LU_C,LU_R)
        ! I made this to use with the Gauss-Seidel method but it turns out I didn't need it.
        ! Might as well keep it...
        use module_common_data , only : p2

        implicit none

        integer, intent(in) :: ncells
        integer, intent(in) :: nnz
        real(p2), dimension(5,5,nnz), intent(in) :: V
        integer, dimension(nnz), intent(in) :: C
        integer, dimension(ncells+1), intent(in) :: R

        real(p2), dimension(:,:,:),allocatable, INTENT(OUT) :: LU_V
        integer, dimension(:),allocatable, INTENT(OUT) :: LU_C
        integer, dimension(ncells+1), INTENT(OUT) :: LU_R

        integer :: nnz_prime
        integer :: i, j

        allocate(LU_V(5,5,nnz-ncells))
        allocate(LU_C(    nnz-ncells))

        nnz_prime = 1

        do i = 1,ncells
            do j = R(i),(R(i+1)-1)
                if (C(j) /= R(i) ) then ! not on the diagonal
                    LU_V(:,:,nnz_prime) = V(:,:,j)
                    LU_C(nnz_prime) = C(j)
                    nnz_prime = nnz_prime + 1
                end if
            end do
            LU_R(i+1) = nnz_prime
        end do
        
    end subroutine LU_extraction

    function insertion_sort_int(x,length)
        ! Performs insertion sort on the list x of integers
        ! While insertion sort is O(n^2) it is generally very efficient for very short lists. The maximum number of elements this
        ! subroutine will experience is currently 4 so a very simple sorting method should be more efficient than a better scaling
        ! but more complex method.
        
        implicit none
        integer, intent(in)                    :: length
        integer, dimension(length), intent(in) :: x
        integer, dimension(length)             :: insertion_sort_int

        integer :: i,j,x_j
        insertion_sort_int = x

        do i = 1,length
            j = i
            inner : do
                if ( (j <= 1) ) then
                    exit inner
                else if (insertion_sort_int(j-1) < insertion_sort_int(j)) then
                    exit inner
                end if
                x_j = insertion_sort_int(j)
                insertion_sort_int(j) = insertion_sort_int(j-1)
                insertion_sort_int(j-1) = x_j
                j = j-1
            end do inner
        end do
    end function

    subroutine insertion_sort_real(x,index)
        use module_common_data , only : p2
        ! This routine uses insertion sort sort the input vector x and returns the output vector index
        ! which is the index needed to sort x.
        ! Insertion sort is O(n^2) which is generally inefficient.  However its simplicity allows it 
        ! to be faster than most O(n log n) methods for small arrays.  Since the incoming array will 
        ! just be the cell neighbors for the coarse mesh n <= 6 unless I get around to adding support
        ! for general polyhedral cells (please no...)
        !
        ! ref: https://en.wikipedia.org/wiki/Insertion_sort

        real(p2), dimension(:), intent(in)  :: x
        integer, dimension(:), intent(out) :: index
        
        real(p2), dimension(:), allocatable :: x_internal
        integer i, j, length, index_j
        real(p2) x_j
        length = size(x)
        allocate(x_internal(length))
        x_internal = x

        index = (/(i,i=1,length)/)
        
        do i = 1,length
            j = i
            inner : do 
                if ( (j <= 1) ) then
                    exit inner
                else if ( x_internal(j-1) < x_internal(j)) then ! it doesn't like evaluating this line with .or. if j = 1 (array bounds)
                                              ! I know this is stupid, but I don't want to debug it right now...
                    exit inner
                end if
                x_j = x_internal(j)
                x_internal(j) = x_internal(j-1)
                x_internal(j-1) = x_j
                index_j = index(j)
                index(j) = index(j-1)
                index(j-1) = index_j
                j = j-1
            end do inner
        end do
        deallocate(x_internal)
    end subroutine insertion_sort_real

end module module_linear_solver
