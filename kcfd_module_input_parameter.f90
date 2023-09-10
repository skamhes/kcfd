module module_input_parameter

  implicit none

  public
  !-------------------------------------------------------------------------
  ! Define 'p2' for double precidion real .

  integer, parameter :: p2 = selected_real_kind(P=15)

  !-------------------------------------------------------------------------
  ! Define input variables and specify default values
  
  ! Project name
    character(80) ::      project_name = "default"    ! project name
    character(80) ::         grid_type = "ugrid"
  
  ! To write a Tecplot data file
    logical       :: generate_tec_file_b = .true.  ! generate_tec_file = T
    logical       :: generate_tec_file_v = .true.  ! generate_tec_file = T

    
  ! Mach number
    real(p2) ::    M_inf = 0.3_p2        ! Freestream Mach number
    real(p2) ::      aoa = 0.0_p2        ! Angle of attack in degrees (x --> z)
    real(p2) :: sideslip = 0.0_p2        ! Sideslip angle in degrees (x --> y)
  
  ! Scheme/solver parameters
    real(p2)               :: CFL                    = 0.5_p2
    logical                :: CFL_ramp               = .false.
    real(p2)               :: CFL_start              = 0.1_p2
    integer                :: CFL_start_iter         = 10
    integer                :: CFL_steps              = 100
    integer                :: solver_max_itr         = 1000
    real(p2)               :: solver_tolerance       = 1.0e-05_p2
    character(80)          :: inviscid_flux          = "roe"
    character(80)          :: inviscid_jac           = "roe"
    character(80)          :: solver_type            = "rk"
    character(80)          :: jacobian_method        = "analytical"
    real(p2), dimension(5) :: eig_limiting_factor    = (/ 0.1, 0.1, 0.1, 0.1, 0.1 /)  !eigenvalue limiting factor
    real(p2), dimension(5) :: variable_ur            = (/ 1, 1, 1, 1, 1 /)  ! Variable under relaxation factors (only used in 
    ! implicit) computations
  
  ! Linear relaxation (preconditioner)
    character(80) :: lrelax_scheme       = "gs"    ! preconditioner scheme type
    integer       :: lrelax_sweeps       = 500     ! preconditioner max relaxation
    real(p2)      :: lrelax_tolerance    = 0.1_p2  ! preconditioner tolerance (reduction)
    integer       :: max_amg_levels       = 5
  
    integer                :: accuracy_order = 1
    logical                ::  use_limiter = .false.
    logical                :: write_data = .false.
    logical                :: import_data = .false.
    logical                :: perturb_initial = .false.
    logical                :: use_amg = .true. ! use algebraic multigrid to accelerate convergence of linear implicit solver

    ! Navier-Stokes Info
    logical :: navier_stokes = .false.
    character(80) :: visc_flux_method = 'alpha'
    real(p2) :: R = 287.058 ! ideal gas constant for air
    real(p2) :: Freestream_Temp = 293.15 ! degK Sea level room temp
    real(p2) :: Pr = 0.72 ! Prandtl number for sea-level air
    real(p2) :: Reynolds = 1.6e6 ! Freestream Reynolds number
    real(p2) :: C_0 = 110.5 ! degK Sutherland's constant (or something like that...) for air
  
  ! End of Default input values
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  ! Group them into "input_parameters":
  
   namelist / input_parameters / &
    project_name         , &
    generate_tec_file_b  , &
    generate_tec_file_v  , &
    M_inf                , &
    aoa                  , &
    sideslip             , &
    inviscid_flux        , &
    inviscid_jac         , &
    solver_type          , &
    CFL                  , &
    solver_tolerance     , &
    solver_max_itr       , &
    eig_limiting_factor  , &
    accuracy_order       , &
    eig_limiting_factor  , &
    write_data           , &
    import_data          , &    
    use_limiter          , &
    perturb_initial      , &
    lrelax_scheme        , &
    lrelax_sweeps        , &
    lrelax_tolerance     , &
    variable_ur          , &
    use_amg              , &
    max_amg_levels       , &
    navier_stokes        , &
    R                    , &
    Reynolds             , &
    Freestream_Temp      , &
    CFL_ramp             , &
    CFL_start            , &
    CFL_start_iter       , &
    CFL_steps            , &
    jacobian_method      , &
    grid_type            , &
    visc_flux_method

    contains    
        subroutine read_nml_input_parameters(namelist_file)
          implicit none
          character(len=*), intent(in) :: namelist_file
          integer :: os

          write(*,*)
          write(*,*) "-------------------------------------------------------"
          write(*,*) " Reading the input file: ", namelist_file ,"...."
          write(*,*)

          open(unit=10,file=trim(namelist_file),form='formatted',status='old',iostat=os)
          read(unit=10,nml=input_parameters)
        
          write(*,*)
          write(*,*) " List of given namelist variables and their values"
          write(*,*)
        
          write(*,nml=input_parameters) ! Print the namelist variables on screen.
          close(10)
        
          write(*,*)
          write(*,*) " End of Reading the input file: input.nml..... "
          write(*,*) "-------------------------------------------------------"
          write(*,*)
        
            

        end subroutine read_nml_input_parameters


end module