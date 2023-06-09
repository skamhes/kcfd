module module_steady_solver
    use module_common_data    , only : p2
    implicit none
    !This module contains the following subroutines:

    public :: steady_solve                 ! solve residual equations (nonlinear system)

    public :: compute_residual_norm        ! compute residual norms
    public :: compute_local_time_step_dtau ! compute pseudo time step
    

    !Nonlinear solvers:

    !public :: explicit_pseudo_time_forward_euler  !Explicit forward Euler
    !public :: explicit_pseudo_time_tvdrk          !Explicit 2-stage TVD RK
    !public :: implicit                            !Implicit defect-correction
    !public :: jfnk                                !Jacobian-Free Newton-Krylov
contains
    subroutine steady_solve
        use module_common_data    , only : p2, half, one, zero, i_iteration, du, lrelax_sweeps_actual, lrelax_roc
        use module_input_parameter, only : solver_type, accuracy_order, inviscid_flux, CFL, solver_max_itr, solver_tolerance, &
                                           import_data
        use module_ccfv_data_soln , only : set_initial_solution, u, w, res, dtau, u2w, load_data_file, &
                                                res_norm, res_norm_initial!, gradw, wsn
        use module_ccfv_data_grid , only : cell, ncells!, face
        use module_ccfv_gradient  , only : construct_vertex_stencil, compute_lsq_coefficients
        use module_ccfv_residual  , only : compute_residual
        implicit none

        !integer                       :: i_iteration = 0
        !real(p2), dimension(5)        :: res_norm, res_norm_initial
        
        !real(p2), dimension(5,ncells) :: local_res
        integer                       :: i, n_residual_evaluation
        integer                       :: L1 = 1
        
        
        i_iteration = 0
        if (import_data) then
            call load_data_file
        else
            call set_initial_solution
        end if

        write(*,*) " ---------------------------------------"
        write(*,*) " Begin Nonlinear iterations"
        write(*,*) 
        write(*,*)
        write(*,*) "    solver_type = ", solver_type
        write(*,*) " accuracy_order = ", accuracy_order
        write(*,*) " inviscid_flux  = ", inviscid_flux
        write(*,*) "            CFL = ", CFL

        call compute_lsq_coefficients

        write(*,*) 
        write(*,*)
        if (trim(solver_type) == "implicit") then
            write(*,*) " Iteration   continuity   x-momemtum   y-momentum   z-momentum    energy       max-res", &
                "    |   proj     reduction"
            allocate( du(5,ncells)) ! allocate du only if it needed
        else
            write(*,*) " Iteration   continuity   x-momemtum   y-momentum   z-momentum    energy       max-res"
        end if
        
        ! Quick initialization
        lrelax_sweeps_actual = 0
        lrelax_roc = 0.0_p2
        n_residual_evaluation = 0


        solver_loop : do while (i_iteration <= solver_max_itr)
            call compute_residual
            !local_res = res
            call compute_residual_norm(res_norm)
            if (i_iteration == 0 .and.(.not.import_data)) then
                res_norm_initial = res_norm
                do i = 1,5
                    if (abs(res_norm_initial(i)) < 1e-014_p2) then
                        res_norm_initial(i) = one ! prevents infinity res/res_norm_init
                    end if
                end do
            end if
            if ( trim(solver_type) == "implicit" ) then

                write(*,'(i10,6es13.3,a,i6,es12.1)') i_iteration, res_norm(:), maxval(res_norm(:)/res_norm_initial(:)), &
                                                  "   | ", lrelax_sweeps_actual, lrelax_roc
            else
                write(*,'(i10,6es13.3)') i_iteration, res_norm(:), maxval(res_norm(:)/res_norm_initial(:))
            end if

            if (maxval(res_norm(:)/res_norm_initial(:)) < solver_tolerance) then
                write(*,*) " Solution is converged!"
                exit solver_loop
            end if
            i_iteration = i_iteration + 1

            call compute_local_time_step_dtau

            ! Update u: u = u + du
            if (trim(solver_type) == "rk") then
                call explicit_pseudo_time_rk
            elseif (trim(solver_type) == "implicit") then
                call implicit
            else
                write(*,*) " Unsopported iteration method: Solver = ", solver_type
            end if


        end do solver_loop
    
    end subroutine steady_solve 

    subroutine explicit_pseudo_time_rk
        use module_common_data    , only : p2, half
        use module_ccfv_data_soln , only : u, w, res, dtau, u2w
        use module_ccfv_data_grid , only : cell, ncells!, face
        use module_ccfv_residual  , only : compute_residual
        implicit none
        real(p2), dimension(5,ncells) :: u0!, local_w,local_u
        integer                       :: i
        
        ! 2 Stage Explicit Runge-kutta
        u0 = u
        do i = 1,ncells
            u(:,i) = u0(:,i) - (dtau(i)/cell(i)%vol) * res(:,i) ! u*
            w(:,i) = u2w(u(:,i))
        end do
        !local_u = u
        !local_w = w
        call compute_residual
        !local_res = res
        ! 2nd stage of Runge-Kutta
        do i = 1,ncells
            u(:,i) = half*(u(:,i) + u0(:,i)) - half * (dtau(i)/cell(i)%vol) * res(:,i) ! u*
            w(:,i) = u2w(u(:,i))
        end do
    end subroutine explicit_pseudo_time_rk

    subroutine implicit
        use module_jacobian        , only : compute_jacobian
        use module_common_data     , only : p2, du
        use module_ccfv_data_grid  , only : ncells
        use module_ccfv_data_soln  , only : u, w, u2w
        use module_linear_solver   , only : linear_relaxation
        implicit none
        integer         :: i
        real(p2)        :: omegan !under-relaxation factor for nonlinear iteration
        ! Compute the residual Jacobian
        call compute_jacobian

        ! Compute du by relaxing the linear system
        call linear_relaxation

        loop_cells : do i = 1,ncells
            !Check negative desnity/pressure and compute a safety factor if needed.
            omegan = safety_factor( u(:,i), du(:,i) )

            !Update the conservative variables
            u(:,i) = u(:,i) + omegan*du(:,i)

            !Update the primitive variables.
            w(:,i) = u2w( u(:,i) )
        end do loop_cells
            
    end subroutine implicit
    !********************************************************************************
    !
    ! This subroutine computes the residual L1 norm (average of the absolute value).
    !
    !********************************************************************************
    subroutine compute_residual_norm(res_norm)

        use module_common_data   , only : p2, zero
        use module_ccfv_data_grid, only : ncells
        use module_ccfv_data_soln, only : res
    
        implicit none
    
        real(p2), dimension(5), intent(out) :: res_norm
    
        !Local variables
        integer                :: i
    
        !Initialize the norm:
        res_norm(:) =  zero
    
        cell_loop : do i = 1, ncells
    
            res_norm(:) = res_norm(:) + abs( res(:,i) ) !L1 norm
    
        end do cell_loop
    
        res_norm(:) = res_norm(:) / real(ncells,p2)   !L1 norm
  
    end subroutine compute_residual_norm
    !********************************************************************************
    subroutine compute_local_time_step_dtau
        use module_common_data      , only : half
        use module_ccfv_data_grid   , only : ncells, cell
        use module_ccfv_data_soln   , only : dtau, wsn
        use module_input_parameter  , only : CFL

        implicit none

        integer :: i

        cell_loop : do i = 1,ncells
            dtau(i) = CFL * cell(i)%vol/( half * wsn(i) )
        end do cell_loop
    end subroutine compute_local_time_step_dtau
!********************************************************************************
! Compute a safety factor (under relaxation for nonlinear update) to make sure
! the updated density and pressure are postive.
!
! This is just a simple exmaple.
! Can you come up with a better and more efficient way to control this?
!
!********************************************************************************
    function safety_factor(u,du)

        use module_common_data          , only : p2
        use module_ccfv_data_soln       , only : u2w
        use module_input_parameter      , only : M_inf
       
        implicit none
       
        real(p2) ::    zero = 0.00_p2
       
        real(p2), dimension(5), intent(in) :: u, du
        real(p2)                           :: safety_factor
        integer                            :: ir = 1, ip = 5
        real(p2), dimension(5)             :: u_updated
        real(p2), dimension(5)             :: w_updated
        real(p2)                           :: p_updated
       
        ! Default safety_factor
    
        safety_factor = 1.0_p2
    
        ! Temporarily update the solution:
    
        u_updated = u + safety_factor*du
        w_updated = u2w(u_updated)
    
        !-----------------------------
        ! Return if both updated density and pressure are positive
    
        if ( w_updated(ir) > zero .and. w_updated(ip) > zero ) then
    
            !Good. Keep safety_factor = 1.0_p2, and return.
    
            return
    
        endif
    
        !-----------------------------
        ! Negative density fix
    
        if ( w_updated(ir) <= zero) then ! meaning du(ir) < zero, reducing the density.
    
            safety_factor = -u(ir)/du(ir) * 0.25_p2 ! to reduce the density only by half.
    
        endif
    
        !-----------------------------
        ! Negative pressure fix
        !
        ! Note: Take into account a case of density < 0 and pressure > 0.
        !       Then, must check if the safety factor computed above for density
        !       will give positive pressure also, and re-compute it if necessary.
    
        if ( w_updated(ip) <= zero .or. w_updated(ir) <= zero) then
    
            !Note: Limiting value of safety_factor is zero, i.e., no update and pressure > 0.
            do
    
            u_updated = u + safety_factor*du
            w_updated = u2w(u_updated)
            p_updated = w_updated(ip)
    
            ! For low-Mach flows, theoretically, pressure = O(Mach^2).
            ! We require the pressure be larger than 1.0e-05*(free stream Mach)^2.
            ! Is there a better estimate for the minimum pressure?
    
            if (p_updated > 1.0e-05_p2*M_inf**2) exit
    
            safety_factor = 0.75_p2*safety_factor ! reduce the factor by 25%, and continue.
    
            end do
    
        endif
    
    
    end function safety_factor
       
end module module_steady_solver