module module_bc_states
    
    implicit none
    
    public :: get_right_state
    public :: get_viscous_right_state
contains
    subroutine get_right_state(xb,yb,zb,ucL,njk,bc_state_type, ucb)
        
        use module_common_data     , only : p2
        use module_ccfv_data_soln  , only : u2w, w2u
        implicit none

        !Input
        real(p2),                   intent(in) :: xb, yb, zb
        real(p2), dimension(5),     intent(in) :: ucL
        real(p2), dimension(3),     intent(in) :: njk
        character(len=*),           intent(in) :: bc_state_type
        
        !output
        real(p2), dimension(5),    intent(out) :: ucb

        !Local variables
        real(p2), dimension(5) :: wL, wb
        real(p2), dimension(5) :: dummy

        wL = u2w(ucL)

        select case(trim(bc_state_type))
            case('freestream')
                call freestream(wb)
            case('slip_wall')
                call slip_wall(wL,njk,wb)
            case('no_slip_wall')
                call no_slip_wall(wL,njk,wb)
            case('outflow_subsonic')
                call back_pressure(wL,wb)
            case default
                write(*,*) "Boundary condition=",trim(bc_state_type),"  not implemented."
                write(*,*) " --- Stop at get_right_state() in kcfd_module_bc_states.f90..."
                stop
        end select

        ucb = w2u(wb)
    end subroutine get_right_state
    
    subroutine get_viscous_right_state(xb,yb,zb,ucL,gradwL,njk,bc_state_type, ucb, gradwb)
        
        use module_common_data     , only : p2, zero
        use module_ccfv_data_soln  , only : u2w, w2u
        implicit none

        !Input
        real(p2),                   intent(in) :: xb, yb, zb
        real(p2), dimension(5),     intent(in) :: ucL
        real(p2), dimension(3,5),   intent(in) :: gradwL
        real(p2), dimension(3),     intent(in) :: njk
        character(len=*),           intent(in) :: bc_state_type
        
        !output
        real(p2), dimension(5),    intent(out) :: ucb
        real(p2), dimension(3,5),  INTENT(OUT) :: gradwb ! This gradient term feels like more of a physical BC.  I may read how to 
                                                         ! do this more "correctly" in the future...

        !Local variables
        real(p2), dimension(5) :: wL, wb
        real(p2), dimension(5) :: dummy

        wL = u2w(ucL)
        gradwb = zero
        select case(trim(bc_state_type))
            case('freestream')
                call freestream(wb)
            case('slip_wall')
                call slip_wall_visc(wL,njk,wb)
                gradwb(:,2:4) = gradwL(:,2:4)
            case('no_slip_wall')
                call no_slip_wall_visc(wL,njk,wb)
                gradwb(:,2:4) = gradwL(:,2:4)
            case('outflow_subsonic')
                call back_pressure_visc(wL,wb)
                gradwb(:,5) = gradwL(:,5)
            case default
                write(*,*) "Boundary condition=",trim(bc_state_type),"  not implemented."
                write(*,*) " --- Stop at get_right_state() in kcfd_module_bc_states.f90..."
                stop
        end select

        ucb = w2u(wb)
    end subroutine get_viscous_right_state

    subroutine freestream(wb)
        use module_common_data   , only : p2
        use module_ccfv_data_soln, only : rho_inf, u_inf, v_inf, w_inf, p_inf
        implicit none

        real(p2),dimension(5), intent(out) :: wb

        wb(1) = rho_inf
        wb(2) = u_inf
        wb(3) = v_inf
        wb(4) = w_inf
        wb(5) = p_inf

    end subroutine freestream

    subroutine back_pressure(wL,wb)
        use module_common_data   , only : p2
        use module_ccfv_data_soln, only : p_inf
        implicit none
        real(p2), dimension(5), intent( in) :: wL
        real(p2),dimension(5), intent(out) :: wb

        wb = wL
        wb(5) = p_inf !<- Just fix the pressure.

    end subroutine back_pressure

    subroutine back_pressure_visc(wL,wb)
        use module_common_data   , only : p2, two
        use module_ccfv_data_soln, only : p_inf
        implicit none
        real(p2), dimension(5), intent( in) :: wL
        real(p2),dimension(5), intent(out) :: wb

        wb = wL
        wb(5) = two * p_inf - wL(5) !<- Just fix the pressure.
                ! = P_L - 2*(P_L-P_inf)
        ! Ensures the arithmetic mean of pressure at the boundary is always equal to p_inf
    end subroutine back_pressure_visc

    subroutine slip_wall(wL,njk,wb)
        use module_common_data     , only : p2
        use module_input_parameter , only : navier_stokes
        use module_ccfv_data_soln  , only : p_inf
        implicit none

        real(p2), dimension(5), intent( in) :: wL
        real(p2), dimension(3), intent( in) :: njk
        real(p2), dimension(5), intent(out) :: wb

        real(p2) :: un
        
        un = wL(2)*njk(1) + wL(3)*njk(2) + wL(4)*njk(3)
        wb = wL
        ! Ensure zero normal velocity on average:
        wb(2) = wL(2) - un*njk(1)
        wb(3) = wL(3) - un*njk(2)
        wb(4) = wL(4) - un*njk(3)

        
    end subroutine slip_wall

    subroutine slip_wall_visc(wL,njk,wb)
        ! The viscous subroutines compute the values at a goast cell with a cell center reflected across the boundary face
        ! this is different to the invisced flux boundary conditions.  For inviscid flux you are only interested in the 
        ! reconstructed values directly on either side of the face so values are given at the face.  For viscous flux you need the
        ! values  at the cell center so we need to calculate things slightly differently...
        use module_common_data     , only : p2, two
        use module_input_parameter , only : navier_stokes
        use module_ccfv_data_soln  , only : p_inf
        implicit none

        real(p2), dimension(5), intent( in) :: wL
        real(p2), dimension(3), intent( in) :: njk
        real(p2), dimension(5), intent(out) :: wb

        real(p2) :: un
        
        un = two * (wL(2)*njk(1) + wL(3)*njk(2) + wL(4)*njk(3))
        wb = wL
        ! Ensure zero normal velocity on average:
        wb(2) = wL(2) - un*njk(1)
        wb(3) = wL(3) - un*njk(2)
        wb(4) = wL(4) - un*njk(3)

        
    end subroutine slip_wall_visc

    subroutine no_slip_wall(wL,njk,wb)
        ! no slip wall with zero heat flux (adiabatic condition)
        use module_common_data     , only : p2, zero
        implicit none

        real(p2), dimension(5), intent( in) :: wL
        real(p2), dimension(3), intent( in) :: njk
        real(p2), dimension(5), intent(out) :: wb

        real(p2) :: un
        
        un = wL(2)*njk(1) + wL(3)*njk(2) + wL(4)*njk(3)
        wb = wL
        
        wb(2) = zero
        wb(3) = zero
        wb(4) = zero
    end subroutine no_slip_wall

    subroutine no_slip_wall_visc(wL,njk,wb)
        ! no slip wall with zero heat flux (adiabatic condition)
        use module_common_data     , only : p2, zero
        implicit none

        real(p2), dimension(5), intent( in) :: wL
        real(p2), dimension(3), intent( in) :: njk
        real(p2), dimension(5), intent(out) :: wb

        ! real(p2) :: un
        
        ! un = wL(2)*njk(1) + wL(3)*njk(2) + wL(4)*njk(3)
        wb = wL
        
        wb(2:4) = -wL(2:4)
    end subroutine no_slip_wall_visc
end module module_bc_states