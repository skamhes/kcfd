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
    real(p2), dimension(:,:)  , pointer :: q     ! low-mach primitive variables at cells/nodes (no rho since flow is ~incompr.).
    real(p2), dimension(:,:,:), pointer :: gradw !gradients of w at cells/nodes.
    real(p2), dimension(:,:,:), pointer :: gradq !gradients of q at cells/nodes.

    real(p2), dimension(:), allocatable :: Temp  ! Temperature (only used for Navier stokes EQ)
    real(p2), dimension(:), allocatable :: mu    ! Viscosity
    real(p2), dimension(:,:), allocatable :: gradT ! gradient of temperature
    real(p2), dimension(:)    , pointer :: dtau  !pseudo time step
    real(p2), dimension(:)    , pointer :: wsn   !maximum eigenvalue at faces
    real(p2), dimension(:,:,:), pointer :: tau ! Viscous stress tensor
    real(p2), dimension(:)    , pointer :: uR2 ! Low Mach Preconditioning scaling factor

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
        use module_input_parameter, only : navier_stokes, low_mach_correction
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

        if (.not.low_mach_correction) allocate( u(nq,ncells) )
        allocate( w(nq,ncells) )
        if (low_mach_correction) allocate( q(nq,ncells) )

        if (navier_stokes .and. .not.low_mach_correction) allocate(    Temp(ncells) )
        if (navier_stokes) allocate(      mu(ncells) )
        ! if (navier_stokes) allocate( tau(3,3,ncells) )
        ! tau = zero ! Initial shear will be zero
        
        ! Initialization
        if (.not. low_mach_correction) u = zero
        w = zero
        if (low_mach_correction) q = zero

        !------------------------------------------------------
        ! Allocate low mach correction array (if used)
        if (low_mach_correction) then
            allocate(uR2(ncells))
        end if

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
        if (low_mach_correction) then
            allocate(gradq(3,nq,ncells))
            gradq = zero
            allocate( gradw(3,nq,ncells) )
        else
            allocate( gradw(3,nq,ncells) )
            gradw = zero
            if (navier_stokes) then
                allocate( gradT(3,ncells) )
                ! Initialization
                gradT = zero
            end if
        end if

        !------------------------------------------------------
        ! Allocate residual array

        allocate( res(nq,ncells) )

        ! Initialize 
        res = zero

    end subroutine construct_ccfv_soln_data

    subroutine set_initial_solution
        use module_common_data    , only : p2, one, pi
        use module_ccfv_data_grid , only : ncells
        use module_input_parameter, only : M_inf, aoa, sideslip, perturb_initial, navier_stokes, R,C_0,Freestream_Temp,Reynolds, &
                                           low_mach_correction
        
        implicit none
        
        integer                :: i
        real(p2), dimension(5) :: w_initial, q_initial
        real(p2)               :: T_initial, mu_initial
        

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


        T_initial = one ! T* = p*(gamma)/rho*. I do Like CFD Eq. 4.14.20
        mu_initial = (M_inf/Reynolds) * ((one + (C_0/Freestream_Temp))/(T_initial + (C_0/Freestream_Temp))) * (one) ** 1.5_p2
        
        ! Scaled Nondimensionalized mu.
        if (perturb_initial) then
            w_initial(iu) = 0.2_p2
            w_initial(iv)= 0.1_p2
            w_initial(iw)= 0.15_p2
        end if            
        
        if (low_mach_correction) then
            q_initial = w2q(w_initial)
            !Set initial solution by the free stream values
            cell_loop_prim : do i = 1, ncells
                !Store the initial solution
                w(:,i) = w_initial
                
                !Compute and store conservative variables.
                q(:,i) = w2q( w_initial )
                if (navier_stokes) then 
                    ! Temp(i) =  T_initial
                    mu(i)   = mu_initial
                end if
            end do cell_loop_prim
        else    
            !Set initial solution by the free stream values
            cell_loop : do i = 1, ncells
                !Store the initial solution
                w(:,i) = w_initial
                
                !Compute and store conservative variables.
                u(:,i) = w2u( w_initial )
                if (navier_stokes) then 
                    Temp(i) =  T_initial
                    mu(i)   = mu_initial
                end if
            end do cell_loop
        end if

    end subroutine set_initial_solution

    subroutine load_data_file
        use module_ccfv_data_grid , only : ncells
        use module_common_data, only : filename_data, p2, one, pi
        use module_input_parameter, only : M_inf, aoa, sideslip, Reynolds, C_0, Freestream_Temp, low_mach_correction
        
        implicit none

        integer :: i, os
        real(p2), dimension(5) :: ui

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
        if (os /= 0) then
            write(*,*) trim(filename_data), " not found. Stop! "
            stop
        end if
        read(1,*) res_norm(1),res_norm(2),res_norm(3),res_norm(4),res_norm(5)
        read(1,*) res_norm_initial(1),res_norm_initial(2), &
                    res_norm_initial(3),res_norm_initial(4),res_norm_initial(5)
        
        
        do i = 1,ncells
            if (low_mach_correction) then
                read(1,*) ui(1),ui(2),ui(3),ui(4),ui(5)
                q(:,i) = u2q(ui)
                mu(i) = (M_inf/Reynolds) * ( (one + ( C_0/Freestream_Temp ) )/(q(5,i) + ( C_0/Freestream_Temp )) ) ** 1.5_p2
            else
                read(1,*) u(1,i),u(2,i),u(3,i),u(4,i),u(5,i)
                w(:,i)      = u2w(u(:,i))
                Temp(i) = w(5,i)*gamma/w(1,i)
                mu(i) = (M_inf/Reynolds) * ( (one + ( C_0/Freestream_Temp ) )/(Temp(i) + ( C_0/Freestream_Temp )) ) ** 1.5_p2
            end if
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
    function w2u(w_in) result(u_out)

        use module_common_data, only : p2, one, half
    
        implicit none
    
        real(p2), dimension(5), intent(in) :: w_in ! input
        real(p2), dimension(5)             :: u_out !output
    
        u_out(1) = w_in(ir)
        u_out(2) = w_in(ir)*w_in(iu)
        u_out(3) = w_in(ir)*w_in(iv)
        u_out(4) = w_in(ir)*w_in(iw)
        u_out(5) = w_in(ip)/(gamma-one)+half*w_in(ir)*(w_in(iu)*w_in(iu)+w_in(iv)*w_in(iv)+w_in(iw)*w_in(iw))
  
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
   function u2w(u_in) result(w_out)
  
    use module_common_data, only : p2, one, half
  
    implicit none
  
    real(p2), dimension(5), intent(in) :: u_in ! input
    real(p2), dimension(5)             :: w_out !output
  
      w_out(ir) = u_in(1)
      w_out(iu) = u_in(2)/u_in(1)
      w_out(iv) = u_in(3)/u_in(1)
      w_out(iw) = u_in(4)/u_in(1)
      w_out(ip) = (gamma-one)*( u_in(5) - half*w_out(1)*(w_out(2)*w_out(2)+w_out(3)*w_out(3)+w_out(4)*w_out(4)) )


  
   end function u2w
  !--------------------------------------------------------------------------------
  
  !********************************************************************************
  ! Compute Q from U
  !
  ! ------------------------------------------------------------------------------
  !  Input:  u = conservative variables (rho, rho*u, rho*v, rho*w, rho*E)
  ! Output:  q =    primitive variables (  p,     u,     v,     w,     T)
  ! ------------------------------------------------------------------------------
  !
  ! Note:    E = p/(gamma-1)/rho + 0.5*(u^2+v^2+w^2)
  !       -> p = (gamma-1)*rho*E-0.5*rho*(u^2+v^2w^2)
  !       -> T = p*gamma/rho    
  ! 
  !********************************************************************************
   function u2q(u_in) result(q_out)
  
    use module_common_data, only : p2, one, half
    
    implicit none
  
    real(p2), dimension(5), intent(in) :: u_in ! input
    real(p2), dimension(5)             :: q_out !output
  
      q_out(iu) = u_in(2)/u_in(1)
      q_out(iv) = u_in(3)/u_in(1)
      q_out(iw) = u_in(4)/u_in(1)
      
      q_out(1) = (gamma-one)*( u_in(5) - half*u_in(1)*( q_out(2)*q_out(2) + q_out(3)*q_out(3) + q_out(4)*q_out(4) ) )
      q_out(5) = q_out(1)*gamma/u_in(1)
  
   end function u2q
  !--------------------------------------------------------------------------------

    !********************************************************************************
    ! Compute U from W
    !
    ! ------------------------------------------------------------------------------
    !  Input:  q =    primitive variables (  p,     u,     v,     w,     T)
    ! Output:  u = conservative variables (rho, rho*u, rho*v, rho*w, rho*E)
    ! ------------------------------------------------------------------------------
    !
    ! Note: rho*E = p/(gamma-1) + rho*0.5*(u^2 + v^2 + w^2)
    !       rho   = p/(gamma*T)
    !********************************************************************************
    function q2u(q_in) result(u_out)

        use module_common_data, only : p2, one, half
        
        implicit none
    
        real(p2), dimension(5), intent(in) :: q_in ! input
        real(p2), dimension(5)             :: u_out !output
    
        u_out(1) = q_in(1)*gamma / q_in(5)
        u_out(2) = u_out(1)*q_in(2)
        u_out(3) = u_out(1)*q_in(3)
        u_out(4) = u_out(1)*q_in(4)
        u_out(5) = q_in(1)/(gamma-one)+half*u_out(1)*(q_in(iu)*q_in(iu)+q_in(iv)*q_in(iv)+q_in(iw)*q_in(iw))
    
    end function q2u

    !********************************************************************************
  ! Compute Q from W
  !
  ! ------------------------------------------------------------------------------
  !  Input:  u = conservative variables (rho, rho*u, rho*v, rho*w,     p)
  ! Output:  q =    primitive variables (  p,     u,     v,     w,     T)
  ! ------------------------------------------------------------------------------
  !
  ! Note:    E = p/(gamma-1)/rho + 0.5*(u^2+v^2+w^2)
  !       -> p = (gamma-1)*rho*E-0.5*rho*(u^2+v^2w^2)
  !       -> T = rho/(gamma-1)    
  ! 
  !********************************************************************************
   function w2q(w_in) result(q_out)
  
    use module_common_data, only : p2, one, half
    
    implicit none
  
    real(p2), dimension(5), intent(in) :: w_in ! input
    real(p2), dimension(5)             :: q_out !output
  
      q_out(iu) = w_in(2)/w_in(1)
      q_out(iv) = w_in(3)/w_in(1)
      q_out(iw) = w_in(4)/w_in(1)
      
      q_out(1) = w_in(5)
      q_out(5) = q_out(1)*gamma/w_in(1)
  
   end function w2q
  !--------------------------------------------------------------------------------

    !********************************************************************************
    ! Compute U from W
    !
    ! ------------------------------------------------------------------------------
    !  Input:  q =    primitive variables (  p,     u,     v,     w,     T)
    ! Output:  u = conservative variables (rho, rho*u, rho*v, rho*w, rho*E)
    ! ------------------------------------------------------------------------------
    !
    ! Note: rho*E = p/(gamma-1) + rho*0.5*(u^2 + v^2 + w^2)
    !       rho   = p/(gamma*T)
    !********************************************************************************
    function q2w(q_in) result(w_out)

        use module_common_data, only : p2, one, half
        
        implicit none
    
        real(p2), dimension(5), intent(in) :: q_in ! input
        real(p2), dimension(5)             :: w_out !output
    
        w_out(1) = q_in(1)*gamma / q_in(5)
        w_out(2) = q_in(2)
        w_out(3) = q_in(3)
        w_out(4) = q_in(4)
        w_out(5) = q_in(1)
    
    end function q2w
end module module_ccfv_data_soln