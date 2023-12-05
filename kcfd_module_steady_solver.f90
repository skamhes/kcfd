module module_steady_solver
    use module_common_data    , only : p2
    implicit none
    !This module contains the following subroutines:

    public :: steady_solve                 ! solve residual equations (nonlinear system)

    public :: compute_residual_norm        ! compute residual norms
    public :: compute_local_time_step_dtau ! compute pseudo time step
    
    private
    real(p2), dimension(5,5) :: var_ur_array

    !Nonlinear solvers:

    !public :: explicit_pseudo_time_forward_euler  !Explicit forward Euler
    !public :: explicit_pseudo_time_tvdrk          !Explicit 2-stage TVD RK
    !public :: implicit                            !Implicit defect-correction
    !public :: jfnk                                !Jacobian-Free Newton-Krylov
contains
    subroutine steady_solve
        use module_common_data    , only : p2, half, one, zero, i_iteration, du, lrelax_sweeps_actual, lrelax_roc
        use module_input_parameter, only : solver_type, accuracy_order, inviscid_flux, CFL, solver_max_itr, solver_tolerance, &
                                           import_data, variable_ur, use_limiter, CFL_ramp, CFL_start, CFL_steps, CFL_start_iter
        use module_ccfv_data_soln , only : set_initial_solution, u, w, res, dtau, u2w, load_data_file, &
                                                res_norm, res_norm_initial!, gradw, wsn
        use module_ccfv_data_grid , only : cell, ncells!, face
        use module_ccfv_gradient  , only : construct_vertex_stencil, compute_lsq_coefficients
        use module_ccfv_residual  , only : compute_residual
        use module_ccfv_limiter   , only : phi
        implicit none

        integer                       :: i, n_residual_evaluation
        integer                       :: L1 = 1

        ! Timing Variables
        real                          :: time, totalTime
        real, dimension(2)            :: values
        integer                       :: minutes, seconds

        ! Stop file
        logical                       :: stop_me
        integer                       :: ierr

        real(p2)                      :: CFL_multiplier, CFL_end, CFL_running_mult
        
        var_ur_array = zero
        do i = 1,5
            var_ur_array(i,i) = variable_ur(i)
        end do

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

        if (CFL_ramp) then
            CFL_end = CFL
            CFL = CFL_start
            CFL_multiplier = (CFL_end/CFL_start)**(one/CFL_steps)
            write(*,*) 'CFL Ramping Enabled'
            write(*,*) 'CFL_start: ', CFL_start
            write(*,*) 'CFL_steps: ', CFL_steps
            write(*,*) 'CFL_mult:  ', CFL_multiplier
            CFL_running_mult = CFL_multiplier
            write(*,*)
            write(*,*)
        end if


        if (import_data) then
            write(*,*)  'Continuing iterations from previous solution...'  ! so I know whether or not I'm loading an old file
            write(*,*)
            write(*,*)
        end if


        if (trim(solver_type) == "implicit") then
            write(*,*) " Iteration   continuity   x-momemtum   y-momentum   z-momentum    energy       max-res", &
                "    |   proj     reduction       time    CFL"
            allocate( du(5,ncells)) ! allocate du only if it needed
        else
            write(*,*) " Iteration   continuity   x-momemtum   y-momentum   z-momentum    energy       max-res    |   time     CFL"
        end if
        
        ! Quick initialization
        lrelax_sweeps_actual = 0
        lrelax_roc = 0.0_p2
        n_residual_evaluation = 0
        if (use_limiter) then
            allocate(phi(ncells))
        end if

        solver_loop : do while (i_iteration <= solver_max_itr)
            call compute_residual
            !local_res = res
            call compute_residual_norm(res_norm)
            call dtime(values,time)
            totalTime = time * real(solver_max_itr-i_iteration) ! total time remaining in seconds
            minutes = floor(totalTime/60.0)
            seconds = mod(int(totalTime),60)
            if (i_iteration == 0 .and.(.not.import_data)) then
                res_norm_initial = res_norm
                minutes = 0
                seconds = 0
                do i = 1,5
                    if (abs(res_norm_initial(i)) < 1e-016_p2) then
                        res_norm_initial(i) = one ! prevents infinity res/res_norm_init
                    end if
                end do
            else if (i_iteration <= 5 .and. (.not. import_data)) then
                
                do i = 1,5
                    if (res_norm_initial(i) == one) then
                        if (abs(res_norm(i)) > 1e-014_p2) then
                            res_norm_initial(i) = res_norm(i) ! prevents infinity res/res_norm_init
                        end if
                    end if
                end do

            end if
            

            if ( trim(solver_type) == "implicit" ) then

                write(*,'(i10,6es13.3,a,i6,es12.1,i10.2,a,i2.2,es13.3)') i_iteration, res_norm(:), & 
                                                  maxval(res_norm(:)/res_norm_initial(:)), &
                                                  "   | ", lrelax_sweeps_actual, lrelax_roc, minutes, ":", seconds, CFL
            else
                write(*,'(i10,6es13.3,i10.2,a,i2.2,es13.3)') i_iteration, res_norm(:), maxval(res_norm(:)/res_norm_initial(:)), &
                                                      minutes, ":", seconds, CFL
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
            elseif (trim(solver_type) == 'explicit') then
                call explicit_pseudo_time_forward_euler
            elseif (trim(solver_type) == "implicit") then
                call implicit
            else
                write(*,*) " Unsopported iteration method: Solver = ", solver_type
            end if

            if (CFL_ramp .and. (i_iteration < CFL_steps + CFL_start_iter) .and. i_iteration > CFL_start_iter) then
                CFL_running_mult = CFL_running_mult * CFL_multiplier
                CFL = CFL_start * CFL_running_mult
            elseif (CFL_ramp .and. (i_iteration == CFL_steps + CFL_start_iter)) then
                CFL = CFL_end
            end if

            ! check for stop file
            ! stop file can be created by typing "echo > kcfdstop" in the working directory
            inquire (file = 'kcfdstop', exist = stop_me)
            if (stop_me) then
                write(*,*) "kcfdstop file found! Stopping iterations!"
                open(10,file = 'kcfdstop',status='old',iostat=ierr)
                if (ierr == 0) then
                    close(10,status ='delete',iostat = ierr) ! delete the file so next time can run again
                    if (ierr == 0) then
                        write(*,*) 'kcfdstop successfully deleted!'
                    end if
                end if
                exit solver_loop
            end if
        end do solver_loop
    
    end subroutine steady_solve 

    subroutine explicit_pseudo_time_rk
        use module_common_data    , only : p2, half, one, zero
        use module_ccfv_data_soln , only : u, w, res, dtau, u2w, gamma, Temp, uR2, u2q, q2u, q, q2w
        use module_ccfv_data_grid , only : cell, ncells!, face
        use module_ccfv_residual  , only : compute_residual
        use module_input_parameter, only : low_mach_correction
        use module_gewp,            only : gewp_solve
        
        implicit none
        real(p2), dimension(5,ncells) :: u0, q0!, local_w,local_u
        integer                       :: i, os
        real(p2) :: H, rho_p, rho_T, theta, rho
        real(p2), dimension(5,5) :: preconditioner, dwdq, pre_inv
        ! 2 Stage Explicit Runge-kutta
        if (low_mach_correction) then
            q0 = q
        else
            u0 = u
        end if

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

                call gewp_solve(preconditioner, 5, pre_inv, os)

                if (os .ne. 0) then
                    write(*,*) 'Error inverting precondition matrix at cell: ', i,' Stop!'
                    stop
                end if
                
                q(:,i) = q0(:,i) - (dtau(i)/cell(i)%vol) * matmul( pre_inv,res(:,i) )
                w(:,i) = q2w(q(:,i))
            else
                u(:,i) = u0(:,i) - (dtau(i)/cell(i)%vol) * res(:,i) ! u*
                w(:,i) = u2w(u(:,i))
            end if
            
        end do

        call compute_residual

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

                call gewp_solve(preconditioner, 5, pre_inv, os)
                if (os .ne. 0) then
                    write(*,*) 'Error inverting precondition matrix at cell: ', i,' Stop!'
                    stop
                end if
                ! q0 = u2q( u0(:,i) ) ! Intermediate solution at cell I in primative values Q
                ! qi = u2q(  u(:,i) ) ! Primative vars 

                q(:,i) = half * (q(:,i) + q0(:,i)) - half * (dtau(i)/cell(i)%vol) * matmul( pre_inv,res(:,i) )
                w(:,i) = q2w(q(:,i))
            else
                u(:,i) = half * (u(:,i) + u0(:,i)) - half * (dtau(i)/cell(i)%vol) * res(:,i) !
                w(:,i) = u2w(u(:,i))
            end if
            
        end do

    end subroutine explicit_pseudo_time_rk

    subroutine explicit_pseudo_time_forward_euler
        use module_common_data    , only : p2, half, one, zero
        use module_ccfv_data_soln , only : u, w, res, dtau, u2w, gamma, Temp, uR2, u2q, q2u,q, q2w
        use module_ccfv_data_grid , only : cell, ncells!, face
        use module_ccfv_residual  , only : compute_residual
        use module_input_parameter, only : low_mach_correction
        use module_gewp,            only : gewp_solve
        implicit none

        real(p2), dimension(5) :: qi, dq, du
        integer i, os
        real(p2) :: H, rho_p, rho_T, theta, Tempi, v2, c2, rho
        real(p2), dimension(5,5) :: preconditioner, dwdq, pre_inv, pletcher
        real(p2) :: a

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
                

                call gewp_solve(preconditioner, 5, pre_inv, os)
                if (os .ne. 0) then
                    write(*,*) 'Error inverting precondition matrix at cell: ', i,' Stop!'
                    stop
                end if
                
                dq = - (dtau(i)/cell(i)%vol) * matmul(pre_inv,res(:,i))
                
                q(:,i) = q(:,i) + dq
                w(:,i) = q2w(q(:,i))
            else
                du = - (dtau(i)/cell(i)%vol) * res(:,i)
                
                u(:,i) = u(:,i) + du
                w(:,i) = u2w(u(:,i))
            end if
            
        end do


    end subroutine explicit_pseudo_time_forward_euler

    subroutine implicit
        use module_jacobian        , only : compute_jacobian
        use module_input_parameter , only : jacobian_method, low_mach_correction
        use module_numerical_jacobian, only:compute_numerical_jacobian
        use module_common_data     , only : p2, du
        use module_ccfv_data_grid  , only : ncells
        use module_ccfv_data_soln  , only : u, w, u2w, u2q, q2u, q, q2w
        use module_linear_solver   , only : linear_relaxation
        implicit none
        integer         :: i
        real(p2)        :: omegan !under-relaxation factor for nonlinear iteration
        
        real(p2), dimension(5)        :: du_local, u_local

        ! Compute the residual Jacobian
        if (trim(jacobian_method) == "analytical") then
            call compute_jacobian
        else if (trim(jacobian_method) == "numerical") then
            call compute_numerical_jacobian
        else
            write(*,*) "unsupported jacobian method. Stop!"
            stop
        end if

        ! Compute du by relaxing the linear system
        call linear_relaxation

        loop_cells : do i = 1,ncells
            if (low_mach_correction) then
                ! du <=> dq
                omegan = safety_factor_primative(q(:,i),du(:,i))
                q(:,i) = q(:,i) + omegan*matmul(var_ur_array,du(:,i))
                ! q = u2q( u(:,i) )
                ! dq = du(:,i) ! in this case q <---> u (sorta... Probably gonna wanna rewrite this at some point...)

                ! q = q + dq
                ! u_local = q2u(q)
                ! du_local = u_local - u(:,i)

                ! !Check negative desnity/pressure and compute a safety factor if needed.
                ! omegan = safety_factor( u(:,i), du_local )

                ! ! Update the conservative variables
                ! u(:,i) = u(:,i) + omegan*matmul(var_ur_array,du_local)
                ! ! u(:,i) = u(:,i) + matmul(var_ur_array,du_local)
                ! ! u(:,i) = q2u(q)
                !Update the primitive variables.
                w(:,i) = q2w( q(:,i) )
            else
                omegan = safety_factor( u(:,i), du(:,i) )
    
                ! Update the conservative variables
                u(:,i) = u(:,i) + omegan*matmul(var_ur_array,du(:,i))
                !Update the primitive variables.
                w(:,i) = u2w( u(:,i) )
            end if

            
            
            

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
        use module_common_data      , only : half, p2
        use module_ccfv_data_grid   , only : ncells, cell
        use module_ccfv_data_soln   , only : dtau, wsn
        use module_input_parameter  , only : CFL

        implicit none

        integer :: i
        ! real(p2), dimension(ncells) :: viscous_dtau

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

    !********************************************************************************
! Compute a safety factor (under relaxation for nonlinear update) to make sure
! the updated density and pressure are postive.
!
! This is just a simple exmaple.
! Can you come up with a better and more efficient way to control this?
!
!********************************************************************************
    function safety_factor_primative(q,dq)

        use module_common_data          , only : p2
        use module_ccfv_data_soln       , only : u2w
        use module_input_parameter      , only : M_inf

       
        implicit none
       
        real(p2) ::    zero = 0.00_p2
       
        real(p2), dimension(5), intent(in) :: q, dq
        real(p2)                           :: safety_factor_primative
        integer                            :: ir = 1, ip = 5
        real(p2), dimension(5)             :: q_updated
        real(p2), dimension(5)             :: w_updated
        real(p2)                           :: p_updated
       
        ! Default safety_factor
    
        safety_factor_primative = 1.0_p2
    
        ! Temporarily update the solution:
    
        q_updated = q + safety_factor_primative*dq
        
        !-----------------------------
        ! Return if both updated density and pressure are positive
    
        if ( q_updated(1) > zero .and. q_updated(5) > zero ) then
    
            !Good. Keep safety_factor = 1.0_p2, and return.
    
            return
    
        endif
    
        !-----------------------------
        ! Negative Temperature fix
    
        if ( q_updated(5) <= zero) then ! meaning du(ir) < zero, reducing the density.
    
            safety_factor_primative = -q(5)/dq(5) * 0.25_p2 ! to reduce the density only by half.
    
        endif
    
        !-----------------------------
        ! Negative pressure fix
        !
        ! Note: Take into account a case of density < 0 and pressure > 0.
        !       Then, must check if the safety factor computed above for density
        !       will give positive pressure also, and re-compute it if necessary.
    
        if ( q_updated(1) <= zero .or. q_updated(5) <= zero) then
    
            !Note: Limiting value of safety_factor is zero, i.e., no update and pressure > 0.
            do
    
            q_updated = q + safety_factor_primative*dq
            p_updated = q_updated(1)
    
            ! For low-Mach flows, theoretically, pressure = O(Mach^2).
            ! We require the pressure be larger than 1.0e-05*(free stream Mach)^2.
            ! Is there a better estimate for the minimum pressure?
    
            if (p_updated > 1.0e-05_p2*M_inf**2) exit
    
            safety_factor_primative = 0.75_p2*safety_factor_primative ! reduce the factor by 25%, and continue.
    
            end do
    
        endif
    
    
    end function safety_factor_primative
       
end module module_steady_solver