module module_common_data
    implicit none
     !This module contains the following subroutines:
    ! - set_filenames
    ! - read_grid
    
    ! Everything is in public and accessible to other modules
    public

    !--------------------------------------------------------------------
    ! Constants
    integer , parameter :: p2 = selected_real_kind(P=15) !Double precision
    real(p2), parameter ::   zero = 0.0_p2
    real(p2), parameter ::   half = 0.5_p2
    real(p2), parameter ::    one = 1.0_p2
    real(p2), parameter ::    two = 2.0_p2
    real(p2), parameter ::  three = 3.0_p2
    real(p2), parameter ::   four = 4.0_p2
    real(p2), parameter ::   five = 5.0_p2
    real(p2), parameter ::  third = 1.0_p2/3.0_p2
    real(p2), parameter :: fourth = 1.0_p2/4.0_p2
    real(p2), parameter ::  sixth = 1.0_p2/6.0_p2
    real(p2), parameter :: my_eps = epsilon(one)  !Machine zero w.r.t. 1.0.
    real(p2), parameter :: pi = 3.141592653589793238_p2
    integer             :: ix = 1, iy = 2, iz = 3

    !--------------------------------------------------------------------
    !--------------------------------------------------------------------
    ! File names

    character(80) :: filename_grid         ! input grid filename (.ugrid)
    character(80) :: filename_bc           ! input bc   filename (.bcmap)
    character(80) :: filename_tecplot_b    ! output tecplot boundary filename (.dat)
    character(80) :: filename_tecplot_v    ! output tecplot volume filename (.dat)
    character(80) :: filename_data         ! Output of U array
    
    !--------------------------------------------------------------------
    ! Data used to read a grid file and a BC file.

     !--------------------------------------------------------------------

    !------------------------------------------
    !>> Node data
    integer                             :: nnodes
    real(p2), dimension(:  )  , pointer :: x, y, z

    !------------------------------------------
    !>> Boundary element data
    integer                              :: nb      !# of boundary segments
    !integer      , dimension(:), pointer :: nbnodes !# of boundary nodes
    !integer      , dimension(:), pointer ::  bnode  !List of boundary nodes
    character(80), dimension(:), pointer :: bc_type !type of boundary condition

    !------------------------------------------
    !>> Element conenctivity data

    !>> Triangular element data
    integer                           :: ntria !# of triangles
    integer , dimension(:,:), pointer ::  tria !List of vertices

    !>> Quadrilateral element data
    integer                           :: nquad !# of quadrilaterals
    integer , dimension(:,:), pointer ::  quad !List of vertices
    
    !>> Tetrahedral element data
    integer                           ::  ntet !# of quadrilaterals
    integer , dimension(:,:), pointer ::   tet !List of vertices

    !>> Tetrahedral element data
    integer                           ::  npyr !# of quadrilaterals
    integer , dimension(:,:), pointer ::   pyr !List of vertices

    !>> Tetrahedral element data
    integer                           ::  nprs !# of quadrilaterals
    integer , dimension(:,:), pointer ::   prs !List of vertices

    !>> Tetrahedral element data
    integer                           ::  nhex !# of quadrilaterals
    integer , dimension(:,:), pointer ::   hex !List of vertices

    !>> Linear solver vars
    integer                           :: lrelax_sweeps_actual
    real(p2)                          :: lrelax_roc
    integer                           :: i_iteration
    real(p2), dimension(:,:), allocatable :: du
    !Note: These variables are available within the entire module, and within
    !      so all subroutines contained below.
    !--------------------------------------------------------------------
    !
    !********************************************************************************
!
! Define file names.
!
! Note: Actual filename character variables are defined in "module_common_data".
!       We access them by the "use" statement and assign the names here.
!
!********************************************************************************
    contains
     subroutine set_filenames

        !To access a user defined character "project_name"
        use module_input_parameter       , only : project_name
    
        implicit none
    
        write(*,*)
        write(*,*) "-------------------------------------------------------"
        write(*,*) " Setting up file names..... "
        write(*,*)
    
        !-----------------------------------------------------------------------
        ! Input grid file (.ugrid):
        ! E.g., filename_grid = "test.grid" if project_name = "test".
    
            filename_grid   = trim(project_name) // '.ugrid'
    
            write(*,'(a28,a28)') "   filename_grid = ", trim(filename_grid)
    
        !-----------------------------------------------------------------------
        ! Input boundary condition file (ASCII file)
        ! E.g., filename_bc = "test.bc" if project_name = "test".
    
            filename_bc   = trim(project_name) // '.bcmap'
    
            write(*,'(a28,a28)') "   filename_bc = ", trim(filename_bc)
    
        !-----------------------------------------------------------------------
        ! Output: Tecplot boundary file (ASCII file)
        ! E.g., filename_tecplot = "test_tec.dat" if project_name = "test".
    
            filename_tecplot_b = trim(project_name) // '_b_tec.dat'
    
            write(*,'(a28,a28)') " filename_tecplot_b = ", trim(filename_tecplot_b)

        !-----------------------------------------------------------------------
        ! Output: Tecplot volume file (ASCII file)
        ! E.g., filename_tecplot = "test_tec.dat" if project_name = "test".
    
            filename_tecplot_v = trim(project_name) // '_v_tec.dat'
    
            write(*,'(a28,a28)') " filename_tecplot_b = ", trim(filename_tecplot_b)


        !-----------------------------------------------------------------------
        ! Output: Tecplot boundary file (ASCII file)
        ! E.g., filename_tecplot = "test_tec.dat" if project_name = "test".
    
            filename_data = trim(project_name) // '.kdat'
    
            write(*,'(a28,a28)') "       filename_data = ", trim(filename_data)

    
        write(*,*)
        write(*,*) " End of Setting up file names..... "
        write(*,*) "-------------------------------------------------------"
        write(*,*)
    
    end subroutine set_filenames
    
    subroutine read_grid

        !********************************************************************************
        !* Read the grid and boundary condition file.
        !* Note: This is a combo of the 2D cell centered and 3D node centered grids
        !        which will be fun...
        !* Note: Currently, the EDU3D-Euler works only for pure tetrahedral grids.
        !*       I think you can modify the code to make it work for other elements.
        !*
        !********************************************************************************
        !*
        !* 1. "datafile_grid_in": .ugrid file name. 
        !*
        !*
        !* 2. "datafile_bcmap_in" is assumed have been written in the following format:
        !*
        !*   -----------------------------------------------------------------------
        !*    write(*,*) nb (# of boundaries)
        !*    do i = 1, nb
        !*     write(*,*) i, bc_type
        !*    end do
        !*   -----------------------------------------------------------------------
        !*
        !*   NOTE: bc_type is the name of the boundary condition, e.g.,
        !*
        !*         1. "freestream"
        !*             Roe flux with freestream condition on the right state.
        !*
        !*         2. "outflow_supersonic"
        !*             Roe flux with the interior state as the right state.
        !*             (equivalent to the interior-extrapolation condition.)
        !*
        !*         3. "slip_wall"
        !*             Right state with zero mass flux through the boundary.
        !*
        !*         4. "outflow_subsonic_p0"
        !*             Fix the back pressure. This should work for subsonic flows in a
        !*             large enough domain.
        !*
        !*         !You can implement more BCs.
        !*
        !*
        !********************************************************************************
        !* Data to be read and stored:
        !*
        !* 1. Some numbers
        !*    nnodes        = number of nodes
        !*    ntria         = number of triangular boundary elements
        !*    nquad         = number of quadrilateral boundary elements
        !*    ntet          = number of tetrahedra
        !*    npyr          = number of pyramids
        !*    nprs          = number of prims
        !*    nhex          = number of hexahedra
        !*
        !* 2. Boundary element data:
        !*    tria = list of vertices of each triangular boundary element
        !*    quad = list of vertices of each quadrilateral boundary element
        !*
        !* 3. Volume element data:
        !*    tet  = list of vertices of each tetrahedron
        !*    pyr  = list of vertices of each pyramid
        !*    prs  = list of vertices of each prism
        !*    hex  = list of vertices of each hexehedron
        !*
        !* 4. x, y, z coordinates of ndoes
        !*    x    = x-coordinate of the nodes
        !*    y    = y-coordinate of the nodes
        !*    z    = z-coordinate of the nodes
        !*
        !* 5. Boundary Data:
        !*    nb      = number of boundary groups
        !*    bc_type = boundary condition name
        !*
        !********************************************************************************
        implicit none

        integer :: os
        integer :: i, ncells, dummy_int
        !integer , dimension(100,8) ::   dummy_debug ! use to debug public variables
        write(*,*)
        write(*,*) "-------------------------------------------------------"
        write(*,*) " Reading : ", trim(filename_grid)
        write(*,*)

        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        ! 1. Read the grid file.

        ! Open the input file.
        open(unit=1, file=filename_grid, status="unknown", iostat=os)

        ! Read: get the size of the grid
        read(1,*) nnodes, ntria, nquad, ntet, npyr, nprs, nhex

        ! Write out the grid data.

        write(*,*) " Total grid numbers:"
        write(*,*) "      Nodes = ", nnodes
        write(*,*) "  Triangles = ", ntria
        write(*,*) "      Quads = ", nquad
        write(*,*) "      Tetra = ", ntet
        write(*,*) "      Hexa  = ", nhex
        write(*,*) "   Pyramids = ", npyr
        write(*,*) "     Prisms = ", nprs
        write(*,*)

        ! Allocate node and element arrays

        if (ntria > 0) allocate(tria(ntria,4))
        if (nquad > 0) allocate(quad(nquad,5))
        if (ntet > 0)  allocate(tet( ntet, 4))
        if (npyr > 0)  allocate(pyr( npyr, 5))
        if (nprs > 0)  allocate(prs( nprs, 6))
        if (nhex > 0)  allocate(hex( nhex, 8))
        
        allocate(x(nnodes),y(nnodes),z(nnodes))

        ! READ: Read the nodal coordinates

        write(*,"(A)",advance="no") " Reading nodes..."

        node_read: do i = 1, nnodes
            read(1,*) x(i), y(i), z(i)
        end do node_read
        write(*,"(A)") "...done"

        write(*,"(A)",advance="no") " Building mesh."
        ! Read element-connectivity information
        ! NOTE: EVERY CELL TYPE MUST BE READ IN A SPECIFIC ORDER

        ! Triangles: assumed that the vertices are ordered counterclockwise
        !
        !         v3
        !         /\
        !        /  \
        !       /    \
        !      /      \
        !     /        \
        !    /__________\
        !   v1           v2

        ! READ: read connectivity info for triangles
        ncells = 1
        if ( ntria > 0) then
            do
                if (ncells > ntria) EXIT
                read(1,*) tria(ncells,1),tria(ncells,2),tria(ncells,3)
                ncells = ncells + 1
            end do
        end if
        write(*,"(A)",advance="no") "." ! Progress bar...

        ! Quads: assumed that the vertices are ordered counterclockwise
        !
        !        v4________v3
        !         /        |
        !        /         |
        !       /          |
        !      /           |
        !     /            |
        !    /_____________|
        !   v1             v2

        ! READ: read connectivity info for quadrilaterals
        ncells = 1
        if ( nquad > 0) then
            do
                if (ncells > nquad) EXIT
                read(1,*) quad(ncells,1),quad(ncells,2),quad(ncells,3),quad(ncells,4)
                ncells = ncells + 1
            end do
        end if
        write(*,"(A)",advance="no") "."
        
        ! Read Tria Boundary Group Number
        ncells = 1
        if ( ntria > 0) then
            do
                if (ncells > ntria) EXIT
                read(1,*) tria(ncells,4)
                ncells = ncells + 1
            end do
        end if
        write(*,"(A)",advance="no") "."
        
        ! Read Quad Boundary Group Number
        ncells = 1
        if ( nquad > 0) then
            do
                if (ncells > nquad) EXIT
                read(1,*) quad(ncells,5)
                ncells = ncells + 1
            end do
        end if
        write(*,"(A)",advance="no") "."

        ! Read connectivity info for tets
        ncells = 1
        if ( ntet > 0) then
            do
                if (ncells > ntet) EXIT
                read(1,*) tet(ncells,1),tet(ncells,2),tet(ncells,3),tet(ncells,4)
                ncells = ncells + 1
            end do
        end if
        write(*,"(A)",advance="no") "."

        ! Read connectivity info for pyrs
        ncells = 1
        if ( npyr > 0) then
            do
                if (ncells > npyr) EXIT
                read(1,*) pyr(ncells,1),pyr(ncells,2),pyr(ncells,3),pyr(ncells,4), &
                pyr(ncells,5)
                ncells = ncells + 1
            end do
        end if
        write(*,"(A)",advance="no") "."

        ! Read connectivity info for prisms
        ncells = 1
        if ( nprs > 0) then
            do
                if (ncells > nprs) EXIT
                read(1,*) prs(ncells,1),prs(ncells,2),prs(ncells,3),prs(ncells,4), &
                prs(ncells,5),prs(ncells,6)
                ncells = ncells + 1
            end do
        end if
        write(*,"(A)",advance="no") "."            

        ! Read connectivity info for hexa
        ncells = 1
        if ( nhex > 0) then
            do
                if (ncells > nhex) EXIT
                read(1,*) hex(ncells,1),hex(ncells,2),hex(ncells,3),hex(ncells,4), &
                hex(ncells,5),hex(ncells,6),hex(ncells,7),hex(ncells,8)
                ncells = ncells + 1
            end do
        end if
        write(*,"(A)") "..done"

        ! Close the file
        close(1)

        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        ! 2. Read the boundary condition data file
        write(*,*)
        write(*,*) "-------------------------------------------------------"
        write(*,*) " Reading the boundary condition file: ", &
            trim(filename_bc)
        write(*,*)

        open(unit=2, file=filename_bc, status="unknown", iostat=os)
        read(2,*) nb
        
        allocate(bc_type(nb))
        do i = 1, nb
            read(2,*) dummy_int, bc_type(i)                               
        end do
        !  Print the data
        do i = 1, nb
            write(*,'(a10,i3,a12,a20)') " boundary", i, "  bc_type = ", trim(bc_type(i))
        end do

        write(*,*)

        close(2)
        ! End of Read the boundary condition data file
        !--------------------------------------------------------------------------------

        write(*,*)
        write(*,*) " Finished Reading : ", trim(filename_grid), " and ",&
            trim(filename_bc)
        write(*,*) "-------------------------------------------------------"
        write(*,*)
    end subroutine read_grid
    



end module module_common_data