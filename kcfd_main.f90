!********************************************************************************
! Karsten-CFD code (kcfd).  If I come up with a funnier name in the future I may
! switch to that... we'll see...
!
! -------------------------------------------------------------------------------
!
! This program does the following
! - reads an unstructured 3D grid (.ugrid, tetra and prism cells currently supported)
! - constructs grid data required by CCFV
! - Constructs solution data required by CCFV
! - sets the initial solution 
! - perform TVD 2nd order Runga-Kutta time stepping (local time-stepping) to solve
!       Euler eqns
! - writes tecplot (.dat) files for the boundary and volume data
!
!
!
! This program is heavily based off of the various edu euler codes developed by 
!       Dr. Katate Masatsuka, which have been very helpful.
!
!
! The author is Karsten Hendrickson
!
! Version 1.0
!
! -------------------------------------------------------------------------------
! Future features I would like to add (in no particular order)
! - Additional solution methods
! - Viscous terms (full Navier-Stokes)
! - RANS Modelling
! - MPI/Prallel computing
! - High order methods
! - Adaptive meshing
! - Multi-grid
! - Probably something else I've forgotten...
!
! -------------------------------------------------------------------------------
! This main program utilizes the following module files.
! (The structure for now will be based on the EDU2D-CCFV code.  I will likely 
!       change it as needed moving forward)
!
!
!  - kcfd_module_input_parameter.f90
!    -> This file contains a module that defines the input data variables, and
!       contains a subroutine that reads an input file named as "input.nml".     
!
!  - kcfd_module_common_data.f90
!    -> This file contains a module that defines common data to be used in the
!       main program and other modules (e.g., basic grid data, file names).
!
!  - kcfd_module_write_files.f90
!    -> This file contains a module that defines subroutiens that write output
!       files: tecplot (.dat).
!
!  - kcfd_module_ccfv_data_grid.f90
!    -> This file contains a module that defines grid data required by CCFV,
!       and contains a subroutine for constructing the grid data.
!
!  - kcfd_module_ccfv_data_soln.f90
!    -> This file contains a module that defines solution data specific to
!       the 2D Euler equations and contains related subroutines.
!
!  - kcfd_module_flux.f90
!    -> This file contains a subroutine that computes a numerical flux.
!
!  - kcfd_module_bc_states.f90
!    -> This file contains a subroutine that computes a BC state for weak BC.
!
!  - kcfd_module_ccfv_residual.f90
!    -> This file contains a subroutine that computes the cccfv residual.
!
!  - kcfd_module_explicit_solver.f90
!    -> This file contains a subroutine for performing a physical/pseudo-time stepping.
!
!  - kcfd_module_gradient.f90
!    -> This file contains a subroutine to compute gradients for second order analysis
!
!
!
!
!********************************************************************************

 program kcfd_main
    ! Access data in module_input_parameter
    ! - read_nml_input_parameters : subroutine to read the input file "input.nml"
    ! -         generate_tec_file : write a Tecplot file if set to be .true.

    use module_input_parameter, only : read_nml_input_parameters, &
                                             generate_tec_file_b, &
                                             generate_tec_file_v, &
                                                    project_name, &
                                                     import_data, &
                                                      write_data
! To access the subroutines: set_filenames and read_grid.

    use module_common_data, only : set_filenames, read_grid

    ! To access the subroutine "construct_ccfv_grid_data" in module_ccfv_data_grid.
    
    use module_ccfv_data_grid , only : construct_ccfv_grid_data
    
    ! To access the subroutines "construct_ccfv_soln_data" and "set_initial_solution"
    ! in module_ccfv_data_soln.
    
    use module_ccfv_data_soln, only : construct_ccfv_soln_data!, w2u, u, w, res                                      
    
    use module_steady_solver
    
    ! To access the subroutines: write_tecplot_file and write_vtk_file
    
    use module_write_files    , only : write_tecplot_file_b, write_tecplot_file_v, write_data_file
                                              
    ! To access the subroutine "explicit_steady_solver" in module_explicit_solver.
    
!    use module_explicit_solver, only : explicit_steady_solver, explicit_unsteady_solver
    
!    use module_ccfv_gradient  , only : compute_lsq_coefficients
    
    implicit none

    write(*,*)
    write(*,*) "----------------------------------------------------------------"
    write(*,*)
    write(*,*) "  This is kcfd Version 1.0"
    write(*,*)
    write(*,*) "----------------------------------------------------------------"
    write(*,*)
    !-------------------------------------------------------------------------------
    ! Read the input parameters

    ! call read_nml_input_parameters("input.nml")
    call read_nml_input_parameters("airfoil_input.nml")
    ! call read_nml_input_parameters("ellipse_input.nml")
    ! call read_nml_input_parameters("1D_shock.nml")
    ! call read_nml_input_parameters("flat_plate.nml")
    ! call read_nml_input_parameters("input_om6.nml")
    ! call read_nml_input_parameters("input_bump.nml")
    !-------------------------------------------------------------------------------
    ! Set file names.

    call set_filenames

    call read_grid

    call construct_ccfv_grid_data

    call construct_ccfv_soln_data

    call steady_solve

    if (write_data) then
        call write_data_file
    end if

    if ( generate_tec_file_b ) then
        call write_tecplot_file_b
    end if
    
    if ( generate_tec_file_v ) then
        call write_tecplot_file_v
    end if




end program kcfd_main