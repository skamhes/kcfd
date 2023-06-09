module module_bc_states
    
    implicit none
    
    public :: get_right_state
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
            case default
                write(*,*) "Boundary condition=",trim(bc_state_type),"  not implemented."
                write(*,*) " --- Stop at get_right_state() in kcfd_module_bc_states.f90..."
                stop
        end select

        ucb = w2u(wb)
    end subroutine get_right_state
    
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

    subroutine slip_wall(wL,njk,wb)
        use module_common_data     , only : p2
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
end module module_bc_states