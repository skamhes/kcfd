module module_ccfv_data_soln
    use module_common_data, only : p2, one, zero ! Grabbing some double precision vars
    implicit none

    !Data that can be used in other modules/subroutines.
    public :: nq, u, w, gradw, res, wsn, dtau
    public :: ir, iu, iv, iw, ip
    public :: gamma
    public :: rho_inf, u_inf, v_inf, w_inf, p_inf
    public :: res_norm, res_norm_initial

    ! Public subroutines:
    public :: construct_ccfv_soln_data
    public :: set_initial_solution
    public :: w2u, u2w
    
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    ! Below are the definition of the data used to implement a cell-centered FV method.
    !------------------------------------------------------------------------------------

    !------------------------------------------
    !>> Solution data
    !------------------------------------------

    integer                             :: nq    !# of eqns/solns (5 for 3D Euler/NS).
    real(p2), dimension(:,:)  , pointer :: u     !conservative variables at cells/nodes.
    real(p2), dimension(:,:)  , pointer :: w     !primitive variables at cells/nodes.
    real(p2), dimension(:,:,:), pointer :: gradw !gradients of w at cells/nodes.

    real(p2), dimension(:)    , pointer :: dtau  !pseudo time step
    real(p2), dimension(:)    , pointer :: wsn   !maximum eigenvalue at faces

    !Just for convenience and clarity in acecssing variables in w or residual components.
    integer                           :: ir = 1  !w(ir) = density
    integer                           :: iu = 2  !w(iu) = x-velocity
    integer                           :: iv = 3  !w(iv) = y-velocity
    integer                           :: iw = 4  !w(iw) = y-velocity
    integer                           :: ip = 5  !w(ip) = pressure

    real(p2), dimension(:,:), pointer :: res     !residual vector
    real(p2), dimension(5)            :: res_norm, res_norm_initial
    !------------------------------------------
    !>> Paramteres
    !------------------------------------------

    real(p2) :: gamma = 1.4_p2 !Ratio of specific heats for air

    !Free stream values: will be set in 'set set_initial_solution' for a given
    !free stream Mach number.
    real(p2) :: rho_inf = one
    real(p2) ::   u_inf = one
    real(p2) ::   v_inf = zero
    real(p2) ::   w_inf = zero
    real(p2) ::   p_inf = one/1.4_p2


    !These data will be allocated for a given grid size, and filled in the
    !following subroutine: construct_ccfv_data.

    !------------------------------------------------------------------------------------
    ! End of the data used to implement a cell-centered FV method.
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------


contains
    subroutine construct_ccfv_soln_data
        use module_common_data   , only : zero   !some p2 values

        use module_ccfv_data_grid, only : ncells !# of cells.
        implicit none
        ! 2D Euler/NS equations = 4 equations:
        !
        ! (1)continuity
        ! (2)x-momentum
        ! (3)y-momentum
        ! (4)z-momentum
        ! (5)energy

        nq = 5

        ! Allocate solution arrays

        allocate( u(nq,ncells) )
        allocate( w(nq,ncells) )

        ! Initialization
        u = zero
        w = zero
        !------------------------------------------------------
        ! Allocate pseudo time step array
        allocate( dtau(ncells) )

        ! Initialization
        dtau = zero

        !------------------------------------------------------
        ! Allocate max wave speed array

        allocate( wsn(ncells) )

        ! Initialization
        wsn = zero

        !------------------------------------------------------
        ! Allocate gradient array

        ! Initialize gradient arrays (just zero here).

        allocate( gradw(3,nq,ncells) )

        ! Initialization
            gradw = zero

        !------------------------------------------------------
        ! Allocate residual array

        allocate( res(nq,ncells) )

        ! Initialize 
        res = zero

    end subroutine construct_ccfv_soln_data

    subroutine set_initial_solution
        use module_common_data    , only : p2, one, pi
        use module_ccfv_data_grid , only : ncells
        use module_input_parameter, only : M_inf, aoa, sideslip, perturb_initial
        implicit none
        
        integer                :: i
        real(p2), dimension(5) :: w_initial
        

        !Set free stream values based on the input Mach number.

        rho_inf = one
        u_inf = M_inf*cos(aoa*pi/180_p2)*cos(sideslip*pi/180_p2)
        v_inf = M_inf*sin(sideslip*pi/180)
        w_inf = M_inf*sin(aoa*pi/180_p2)*cos(sideslip*pi/180_p2)
        p_inf = one/gamma
        
        w_initial(ir) = rho_inf
        w_initial(iu) =   u_inf
        w_initial(iv) =   v_inf
        w_initial(iw) =   w_inf
        w_initial(ip) =   p_inf
        if (perturb_initial) then
            w_initial(iu) = 0.4
            w_initial(iv)= 0.1
            !w_inf = 0.2
        end if                
        !Set initial solution by the free stream values
        cell_loop : do i = 1, ncells
            !Store the initial solution
            w(:,i) = w_initial
            
            !Compute and store conservative variables.
            u(:,i) = w2u( w_initial )
            
        end do cell_loop

    end subroutine set_initial_solution

    subroutine load_data_file
        use module_ccfv_data_grid , only : ncells
        use module_common_data, only : filename_data, p2, one, pi
        use module_input_parameter, only : M_inf, aoa, sideslip
        
        implicit none

        integer :: i, os

        rho_inf = one
        u_inf = M_inf*cos(aoa*pi/180_p2)*cos(sideslip*pi/180_p2)
        v_inf = M_inf*sin(sideslip*pi/180)
        w_inf = M_inf*sin(aoa*pi/180_p2)*cos(sideslip*pi/180_p2)
        p_inf = one/gamma


        write(*,*)
        write(*,*) "-------------------------------------------------------"
        write(*,*) " Reading : ", trim(filename_data)
        write(*,*)

        open(unit=1, file=filename_data, status="unknown", iostat=os)
        read(1,*) res_norm(1),res_norm(2),res_norm(3),res_norm(4),res_norm(5)
        read(1,*) res_norm_initial(1),res_norm_initial(2), &
                    res_norm_initial(3),res_norm_initial(4),res_norm_initial(5)
        
        
        do i = 1,ncells
            read(1,*) u(1,i),u(2,i),u(3,i),u(4,i),u(5,i)
        end do
        close(1)
        
    end subroutine 
    !********************************************************************************
! Compute U from W
!
! ------------------------------------------------------------------------------
!  Input:  w =    primitive variables (rho,     u,     v,     w,     p)
! Output:  u = conservative variables (rho, rho*u, rho*v, rho*w, rho*E)
! ------------------------------------------------------------------------------
!
! Note: E = p/(gamma-1)/rho + 0.5*(u^2+v^2)
!
!********************************************************************************
 function w2u(w) result(u)

    use module_common_data, only : p2, one, half
  
    implicit none
  
    real(p2), dimension(5), intent(in) :: w ! input
    real(p2), dimension(5)             :: u !output
  
     u(1) = w(ir)
     u(2) = w(ir)*w(iu)
     u(3) = w(ir)*w(iv)
     u(4) = w(ir)*w(iw)
     u(5) = w(ip)/(gamma-one)+half*w(ir)*(w(iu)*w(iu)+w(iv)*w(iv)+w(iw)*w(iw))
  
   end function w2u
  !------------------------------------------------------------------------------
  
  !********************************************************************************
  ! Compute U from W
  !
  ! ------------------------------------------------------------------------------
  !  Input:  u = conservative variables (rho, rho*u, rho*v, rho*w, rho*E)
  ! Output:  w =    primitive variables (rho,     u,     v,     w,     p)
  ! ------------------------------------------------------------------------------
  !
  ! Note:    E = p/(gamma-1)/rho + 0.5*(u^2+v^2+w^2)
  !       -> p = (gamma-1)*rho*E-0.5*rho*(u^2+v^2w^2)
  ! 
  !********************************************************************************
   function u2w(u) result(w)
  
    use module_common_data, only : p2, one, half
  
    implicit none
  
    real(p2), dimension(5), intent(in) :: u ! input
    real(p2), dimension(5)             :: w !output
  
      w(ir) = u(1)
      w(iu) = u(2)/u(1)
      w(iv) = u(3)/u(1)
      w(iw) = u(4)/u(1)
      w(ip) = (gamma-one)*( u(5) - half*w(1)*(w(2)*w(2)+w(3)*w(3)+w(4)*w(4)) )
  
   end function u2w
  !--------------------------------------------------------------------------------
  
end module module_ccfv_data_soln