module module_linear_solver

    implicit none
    
    public :: linear_relaxation

contains

    subroutine linear_relaxation
        use module_common_data          , only : p2, lrelax_roc, lrelax_sweeps_actual, one, i_iteration, du, &
                                                 zero, my_eps
        use module_input_parameter      , only : solver_type, lrelax_sweeps, lrelax_tolerance, lrelax_scheme
        use module_ccfv_data_soln       , only : res
        use module_ccfv_data_grid       , only : ncells, cell, jac

        implicit none
        
        real(p2), dimension(5)      :: b
        integer                     :: i, ii, k, os

        ! Linear system residual
        real(p2), dimension(5)      :: linear_res

        ! Residual norms(L1,L2,Linf)
        real(p2), dimension(5,3)    :: linear_res_norm, linear_res_norm_init, linear_res_norm_prev

        ! Rate of convergence
        real(p2) :: roc, roc_previous

        ! Under-relaxation parameter
        real(p2) :: omega_lrelax

        ! Initialize some variables
        lrelax_sweeps_actual = 0 ! to count the actual number of relaxations
        omega_lrelax         = one

        !Open the output file.
        open(unit=1000, file='linear_solve.1000', status="unknown", iostat=os)


        if    ( trim(solver_type) == "implicit" ) then

            write(1000,*) "-------------------------------------------------------------------------------"
            write(1000,*) " Nonlinear Iteration = ", i_iteration
            write(1000,*) "-------------------------------------------------------------------------------"
           
        elseif( trim(solver_type) == "jfnk" ) then
           
        !    write(1000,*) "-------------------------------------------------------------------------------"
        !    write(1000,*) " Nonlinear Iteration = ", i_iteration, " GCR projection = ", gcr_proj_actual
        !    write(1000,*) "-------------------------------------------------------------------------------"
           
        endif
        write(1000,*)
        write(1000,*) " lrelax_scheme = ", trim(lrelax_scheme)
        write(1000,*)

        !--------------------------------------------------------------------------------
        ! 1. Initialize the correction
        du = zero

        !--------------------------------------------------------------------------------
        ! 2. Linear Relaxation (Sweep)
        relax : do ii = 1,lrelax_sweeps
            linear_res_norm(:,1) = zero
            !---------------------------------------------------------
            ! Sequential Gauss-Seidel Relaxation(sweep)
            if (trim(lrelax_scheme) == 'gs') then
                gs_loop : do i = 1,ncells
                    ! Form the right hand side of GS: [ sum( off_diagonal_block*du ) - residual ]
                    b = -res(:,i)
                    gs_nghbr_loop : do k = 1,cell(i)%nnghbrs
                        ! Add RHS from off diagonal terms and du (du = zero to start and will be updated as we go)
                        b = b - matmul(jac(i)%off_diag(:,:,k), du(:,cell(i)%nghbr(k)))
                    end do gs_nghbr_loop 
                    ! Update du by the GS relaxation:
                    !
                    ! e.g., for 3 nghbrs, perform the relaxation in the form:
                    !
                    !                     diagonal block        sum of off-diagonal block contributions
                    !       dUj = omega*{ [V/dtj+dR/dUj]^{-1}*(-[dRj/dU1]*dU1 -[dRj/dU2]*dU2 -[dRj/dU3]*dU3 -Res_j) - dUj }
                    linear_res = matmul(jac(i)%diag_inv(:,:), b) - du(:,i)
                    du(:,i) = du(:,i) + omega_lrelax * linear_res
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
                write(1000,'(a,i10,a,es12.5)') "  after relax ", ii, &
                            " max(L1 norm) = ", maxval(linear_res_norm(:,1))

                if ( maxval(linear_res_norm(:,1)) < my_eps ) then
                    write(1000,*) " Machine zero res reached. Exit GS relaxation. Total sweeps = ", ii
                    write(1000,*) " tolerance_linear = ", lrelax_tolerance
                    lrelax_sweeps_actual = ii
                    exit relax
                endif
            ! After the second Relaxation
            else 
                roc = maxval(linear_res_norm(1:5,1)/linear_res_norm_init(1:5,1))

                !-----------------------------------------------------------------
                ! Print the current linear residual norm in the file 'fort.1000'.
                if (roc < one) then
                    write(1000,'(a,i10,a,es12.5,2x,f6.3,2x,f8.3)')   "  after relax ", ii,          &
                                " max(L1 norm) = ", maxval(linear_res_norm(1:5,1)), omega_lrelax, roc
                else

                    write(1000,'(a,i10,a,es12.5,2x,f6.3,2x,f8.3,a)') "  after relax ", ii,          &
                            " max(L1 norm) = ", maxval(linear_res_norm(1:5,1)), omega_lrelax, roc," <- diverge"
                endif
                !-----------------------------------------------------------------
                ! Tolerance met: Exit

                if (roc < lrelax_tolerance) then
                    write(1000,*)
                    write(1000,*) " Tolerance met. Exit GS relaxation. Total sweeps = ", ii
                    write(1000,*) " tolerance_linear = ", lrelax_tolerance
                    lrelax_sweeps_actual = ii
                    exit relax
                !-----------------------------------------------------------------
                ! If tolerance is NOT met.
                else

                    !---------------------------------------------
                    ! Stop if the linear residual is too small.
             
                    if ( maxval(linear_res_norm(1:5,1)) < my_eps ) then
                        write(1000,*)
                        write(1000,*) " Residuals too small. Exit GS relaxation. Total sweeps = ", ii
                        write(1000,*) " maxval(linear_res_norm(1:5,1)) = ", maxval(linear_res_norm(1:5,1))
                        lrelax_sweeps_actual = ii
                        exit relax
                    endif
            
                    !-------------------------------------------------
                    ! Stop if the maximum number of sweeps is reached.
             
                    if (ii == lrelax_sweeps) then
                        write(1000,*)
                        write(1000,*) " Tolerance not met... sweeps = ", lrelax_sweeps
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
            
        end do relax
        ! End of 3. Linear Relaxation (Sweep)
        !--------------------------------------------------------------------------------

        lrelax_roc = roc
        write(1000,*)
    
        return

        close(1000)


    end subroutine linear_relaxation
end module module_linear_solver
