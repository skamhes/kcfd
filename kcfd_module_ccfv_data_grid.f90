module module_ccfv_data_grid
    use module_common_data, only : p2 ! double precision

    implicit none
    !'public' means these data can be accessed in other program, subroutine/function,
    !and modules. Those not explcitly stated are private to this module, and cannot
    !be accessed by other modules and subroutines/functions.

    !To make the data accesible by other modules/subroutines.
    public :: cell, ncells                                  ! Cell data
    public :: nfaces, face, face_nrml, face_nrml_mag        ! Face data
    public :: bound !nbound                                 ! Boundary data
    public :: heffn, heffc, heffv, heffv_min, heffv_max     ! Effective mesh spacing

    ! Make this subroutine accesible to other modules
    public :: construct_ccfv_grid_data ! Constructs the cell-centered FV data

    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    ! Below are the definition of the data used to implement a cell-centered FV method.
    !------------------------------------------------------------------------------------

    !------------------------------------------------
    !>> Cell data: 1D array including trias + quads
    !------------------------------------------------

    integer :: ncells ! number of cells

    ! custom data type for cell-centered method
    type cc_data_type
        integer                         ::    nvtx ! Number of vertices
        integer, dimension(:), pointer  ::     vtx ! list of vertices
        real(p2)                        :: xc, yc, zc ! cell centroid coordinates
        real(p2)                        ::        vol ! Cell volume
        integer                         ::    nnghbrs ! number of neighbors
        integer, dimension(:),  pointer ::      nghbr ! list of face neighbors
        integer, dimension(:,:),pointer ::    commonv ! list of vertices shared with nghbr
        integer, dimension(:),  pointer ::   ncommonv ! number of common vertices with neighbor
        integer                         ::   cnnghbrs ! number of neighbor groups (only used during AMG)
        integer, dimension(:),  pointer ::     cnghbr ! list of neighbor groups (only used during AMG)
    end type cc_data_type

    ! Cell data array in the custom data type.
    type(cc_data_type), dimension(:), pointer :: cell

    !------------------------------------------
    !>> Face data: interior faces
    !------------------------------------------
    integer                             :: nfaces        ! # of interior faces
    !integer                             :: nfacesides    ! # of sides for each interior face
    integer , dimension(:,:), pointer   :: face          ! interior face data
    real(p2), dimension(:,:), pointer   :: face_nrml     ! face normal (unit vector)
    real(p2), dimension(:)  , pointer   :: face_nrml_mag ! face normal magnitude
    real(p2), dimension(:,:), pointer   :: face_centroid
    
    type bgrid_type
        integer                                 :: nbfaces          ! number of faces for boundary ib
        integer, dimension(:,:), pointer        :: bfaces           ! number of faces
        real(p2), dimension(:,:),  pointer      :: bface_nrml       ! normal vector of boundary face (unit vector)
        real(p2), dimension(:),    pointer      :: bface_nrml_mag   ! normal vector magnitude
        integer, dimension(:), pointer          :: bcell            ! adjacent cell of boundary face
        real(p2), dimension(:,:),  pointer      :: bface_center      ! face center coordinates
    end type bgrid_type

    ! Boundary data array using the custom data type
    type(bgrid_type), dimension(:), pointer :: bound
    type(bgrid_type), dimension(:), pointer :: bound_export
    !------------------------------------------------
    !>> Various effective mesh spacings.
    !------------------------------------------------

    real(p2) ::  heffn     !Effective mesh spacing based on number of nodes.
    real(p2) ::  heffc     !Effective mesh spacing based on number of cells.
    real(p2) ::  heffv     !Effective mesh spacing based on sqrt(volume).
    real(p2) ::  heffv_min !Min effective mesh spacing based on sqrt(volume).
    real(p2) ::  heffv_max !Max effective mesh spacing based on sqrt(volume).


    !These data will be allocated for a given grid size, and filled in the
    !following subroutine: construct_ccfv_data.

    ! kth value numbers for implicit solver
    ! kth_nghbr_of_1 = c2 of a face is the kth neighbor of c1
    ! kth_nghbr_of_2 = c1 of a face is the kth neighbor of c2
    integer, dimension(:), allocatable :: kth_nghbr_of_1
    integer, dimension(:), allocatable :: kth_nghbr_of_2

    ! Jacobian type and var
    type jacobian_type
        real(p2), dimension(5,5)                :: diag     ! diagonal blocks of Jacobian matrix
        real(p2), dimension(:,:,:), allocatable :: off_diag ! off-diagonal blocks
        real(p2), dimension(5,5)                :: diag_inv ! inverse of diagonal blocks
        real(p2), dimension(5)                  :: RHS      ! Right hand side (b) of the linear system
    end type jacobian_type

    type(jacobian_type), dimension(:), allocatable :: jac ! jacobian array
    type(jacobian_type), dimension(16)             :: test ! jacobian array
    !------------------------------------------------------------------------------------
    ! End of the data used to implement a cell-centered FV method.
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------

contains
    subroutine construct_ccfv_grid_data
        !Input: data read from a given grid, and some parameters.
        use module_common_data, only : nnodes                , & !# of nodes
        x, y, z               , & !nodal coords
        ntria, tria           , & !# of trias, and list
        nquad, quad           , & !# of quads, and list
        ntet, tet             , & !# of tets, and list
        npyr                  , & !# of pyramids, and list
        nprs, prs             , & !# of prisms, and list
        nhex                  , & !# of hex, and list
        hex, & !, pyr            , & ! currently unsupported
        nb                    , & !# of boundaries
        ix, iy, iz            , & !ix=1, iy=2, iz=3
        p2                    , & !Double precision
        zero, one, three    !some p2 values        
        
        implicit none

        !------------------------------------------------------------------
        ! Below are some local (temporary) variables used to construct the CCFV data.
        !To store the list of cells around each node.
        type node_type
            integer                        :: nc
            integer, dimension(:), pointer :: c
        end type node_type
        ! Array of custom node-type data
        type(node_type), dimension(:), pointer :: node
        ! Debugging vars
        

        ! Other local variables
        integer, dimension(6)       :: temp
        integer, dimension(6,4)     :: tempv
        integer, dimension(6)       :: tempc
        integer                     ::   i,   j,   k,  kk,  ck
        integer                     ::  v1,  v2,  v3,  v4,  v5,  v6, v7, v8 ! these may have been shifted to specific subroutines
        integer                     :: vk1, vk2, vk3, vk4, vk5, vk6, vk7, vk8, vk 
        integer                     ::  c1,  c2                             ! left and right cells for vol calc
        real(p2)                    ::  x1,  x2,  x3,  x4  
        real(p2)                    ::  y1,  y2,  y3,  y4
        real(p2)                    ::  z1,  z2,  z3,  z4                ! these may have been shifted to specific subroutines
        
        integer                     :: cellCounter                  ! Index variables
        integer                     :: pnvtx                        ! # of nodes for parent cell
        !Variables used in verification
        integer                           :: cell1, cell2
        !real(p2)                          :: vol_domain, vol_domain_cells, xm, ym, volk
        real(p2), dimension(3)            :: sum_bnormal
        real(p2), dimension(:,:), pointer :: sum_face_normal
        real(p2)                          :: machine0_1      !to store machine zero relative to 1
        real(p2)                          :: sum_machine0_1  !to store sum of machine zero
        real(p2)                          :: ref_mag         !to store a reference magnitude
        integer                           :: dummyScalar = 0 
        integer, dimension(5)             :: korder = (/1, 4, 1, 2, 3/)     ! order of prism assoc vertex analysis
        real(p2),dimension(3)             :: x_f  ! for the volume calculation
        real(p2)                          :: dot_result1, dot_result2
        logical                           :: inside ! for checking cell center location
        
        !debugging vars
        !type(cc_data_type), dimension(:), pointer :: local_cell
        !integer , dimension(:,:), pointer   :: local_face
        !real(p2), dimension(:,:), pointer   :: local_face_nrml     ! face normal (unit vector)
        !real(p2), dimension(:)  , pointer   :: local_face_nrml_mag ! face normal magnitude
        !real(p2), dimension(:,:), pointer   :: local_face_centroid
        !type(bgrid_type), dimension(:), pointer :: local_bound
        write(*,*)
        write(*,*) "-------------------------------------------------------"
        write(*,*) "-------------------------------------------------------"
        write(*,*) "-------------------------------------------------------"
        write(*,*) " Construct ccfv grid data and verifty them." 
        write(*,*)

        !------------------------------------------------------------------------------------
        !------------------------------------------------------------------------------------
        ! Construct a single custum-data array for cells, which stores cell information
        ! in a single array, where the first ntria elements are for triangles, and
        ! the last nquad elements are for quadrilaterals.
        !------------------------------------------------------------------------------------
        !------------------------------------------------------------------------------------

        !---------------------------------------------------------------
        !---------------------------------------------------------------
        ! (1) Create a single cell array: cell(:)
        !
        !    Note: We already have triangle and quad lists read from
        !          .grid file: tria(:,3) and quad(:,4). We will store
        !          them in a single custom-data array:
        !
        !               cell(i), i=1,ncells.
        !---------------------------------------------------------------

        write(*,*) " Creating cell array....."
        
        ! Count number of cells
        ncells = ntet + nhex + nprs + npyr
        ! check that there are only supported cells (currently tet and prs)
        if (ncells /= ntet + nprs +  nhex) then
            write(*,*) " Unsupported Grid Elements Present. Stop"
            stop
        end if

        ! allocate the single cell array
        allocate( cell(ncells))
        allocate( jac( ncells))
        cellCounter = 1
        
        tet_loop : do i = 1, ntet
            cell(cellCounter)%nvtx = 4
            allocate(cell(cellCounter)%vtx(4))  ! allocate array of size 4
            cell(cellCounter)%vtx(:) = tet(i,:) ! I should switch the rank major of this but I don't want to right now...
            cellCounter = cellCounter + 1
        end do tet_loop

        prs_loop : do i = 1, nprs
            cell(cellCounter)%nvtx = 6
            allocate(cell(cellCounter)%vtx(6))  ! allocate array of size 6
            cell(cellCounter)%vtx(:) = prs(i,:)
            cellCounter = cellCounter + 1
        end do prs_loop

        hex_loop : do i = 1,nhex
            cell(cellCounter)%nvtx = 8
            allocate(cell(cellCounter)%vtx(8)) ! allocate array of size 8
            cell(cellCounter)%vtx(:) = hex(i,:)
            cellCounter = cellCounter + 1
        end do hex_loop
        !Now we can loop over all cells by cell(i), i=1,ncells, instead
        !of loop over triangles and then quads. Good.

        !---------------------------------------------------------------
        !---------------------------------------------------------------
        ! (2) Compute and store the centroid coordiantes.
        !---------------------------------------------------------------
        write(*,*) " Computing centroid coordinates....."
        ! Centroid is calculated as the mean of each node
        cell_loop : do i = 1,ncells 
            if (cell(i)%nvtx == 4) then ! Tet
                v1 = cell(i)%vtx(1)
                v2 = cell(i)%vtx(2)
                v3 = cell(i)%vtx(3)
                v4 = cell(i)%vtx(4)
                cell(i)%xc = ( x(v1) + x(v2) + x(v3) + x(v4) ) / 4.0_p2
                cell(i)%yc = ( y(v1) + y(v2) + y(v3) + y(v4) ) / 4.0_p2
                cell(i)%zc = ( z(v1) + z(v2) + z(v3) + z(v4) ) / 4.0_p2
            elseif (cell(i)%nvtx == 6) then ! prism cell
                v1 = cell(i)%vtx(1)
                v2 = cell(i)%vtx(2)
                v3 = cell(i)%vtx(3)
                v4 = cell(i)%vtx(4)
                v5 = cell(i)%vtx(5)
                v6 = cell(i)%vtx(6)
                cell(i)%xc = ( x(v1) + x(v2) + x(v3) + x(v4) + x(v5) + x(v6) ) / 6.0_p2
                cell(i)%yc = ( y(v1) + y(v2) + y(v3) + y(v4) + y(v5) + y(v6) ) / 6.0_p2
                cell(i)%zc = ( z(v1) + z(v2) + z(v3) + z(v4) + z(v5) + z(v6) ) / 6.0_p2
            elseif (cell(i)%nvtx == 8) then ! hex cell
                v1 = cell(i)%vtx(1)
                v2 = cell(i)%vtx(2)
                v3 = cell(i)%vtx(3)
                v4 = cell(i)%vtx(4)
                v5 = cell(i)%vtx(5)
                v6 = cell(i)%vtx(6)
                v7 = cell(i)%vtx(7)
                v8 = cell(i)%vtx(8)
                cell(i)%xc = ( x(v1) + x(v2) + x(v3) + x(v4) + x(v5) + x(v6) + x(v7) + x(v8) ) / 8.0_p2
                cell(i)%yc = ( y(v1) + y(v2) + y(v3) + y(v4) + y(v5) + y(v6) + y(v7) + y(v8) ) / 8.0_p2
                cell(i)%zc = ( z(v1) + z(v2) + z(v3) + z(v4) + z(v5) + z(v6) + z(v7) + z(v8) ) / 8.0_p2
            else
                write(*,*) " Something is wrong. Invalid cell(i)%nvtx = ", cell(i)%nvtx, &
                        " at cell #: ", i
            end if
        end do cell_loop
        !---------------------------------------------------------------
        !---------------------------------------------------------------
        ! (4) Construct face-neighbor lists
        !---------------------------------------------------------------
        !-------------------------------------------------------
        ! First, create node-to-cell lists for convenience.
        !
        ! Example: node i has 4 cells around it,
        !
        !        o-------o-------------o
        !       /        |             |
        !      /    23   |      41     |
        !     o----------o-------------o
        !      \        i \            |
        !       \   101    \     13    |
        !        \          \          | 
        !         o----------o---------o
        !
        !          and so we'll have
        !            node(i)%nc     = 4
        !            node(i)%c(1:4) = [13,23,41,101] !<- No particular order.
        !
        write(*,*) " Creating node dependency lists..."
        
        allocate(node(nnodes))
        ! Initialize the # of cells around a node.
        do i = 1, nnodes
            node(i)%nc  = 0
        end do
        ! Count the number of cells around a node
        do i = 1, ncells
            do k = 1, cell(i)%nvtx
                vk = cell(i)%vtx(k) ! get node(k) from cell(i)
                node(vk)%nc = node(vk)%nc + 1 ! increase adjacent cell count by one for node(k)
            end do
        end do
        ! We can now allocate the node-to-cell arrays
        do i = 1, nnodes
            allocate( node(i)%c(node(i)%nc))
        end do
        ! Re-initialize the # of cells around the node
        do i = 1, nnodes
            node(i)%nc  = 0
        end do
        ! Fill the node-to-cell arrays
        do i = 1, ncells
            do k = 1, cell(i)%nvtx
                vk = cell(i)%vtx(k) ! node(k) from cell(i)
                node(vk)%nc = node(vk)%nc + 1 ! increase adjacent cell count nc
                node(vk)%c(node(vk)%nc) = i ! for node vk, the nc_th cell is equal to i
            end do
        end do

        !-------------------------------------------------------
        !-------------------------------------------------------
        ! Now create the face-neighbor list for each cell.
        !
        ! Example: The cell i below has 4 face neighbors:
        !          the cells, 3, 4, 6, 9,
        !
        !          / \
        !         / 9 \
        !  _____ /_____\
        !  \  4 /      /\
        !   \  /   i  /  \
        !    \/______/__3_\
        !     \     /
        !      \ 6 /
        !       \ /
        !
        ! and so we'll have
        !
        !         cell(i)%nnghbrs    = 4
        !         cell(i)%nghbr(1:4) = [3,6,9,4] !<- No particular order here.
        !
        ! and for each face, k, we'll also store the nodes(vertices) that
        ! defines the face:
        !
        !         cell(i)%commonv(k,1) !right node seen from i to k
        !         cell(i)%commonv(k,2) ! left node seen from i to k
        !
        !        _______2 (left node)
        !       /      /\
        !      /   i  /k \
        !     /______/____\
        !           1 (right node)
        !  	      c
        ! 	     /\
        ! 	    /| \
        ! 	   / |  \
        ! 	  /  *d  \
        ! 	 / *   *  \
        ! 	/__________\b
        !   a
        ! 
        ! Front: 1 ==> a-b-c
        ! Right: 2 ==> b-c-d
        ! Left:  3 ==> c-d-a
        ! Base:  4 ==> d-a-b

        write(*,*) " Constructing face-neighbor lists....."
        cell_assoc_loop : do i = 1,ncells
            pnvtx = cell(i)%nvtx
            cell(i)%nnghbrs = 0
            ! call associate_face(i)
            ! check for supported cells (currently tet and prs)
            if (.not.(pnvtx == 4 .or. pnvtx == 6 .or. pnvtx == 8)) then
                write(*,*) " Unsopported cell shape"
                stop
            end if

            cell_type : if (pnvtx == 4) then ! parent cell is a tet
                face_k_tet_loop : do k = 1, (cell(i)%nvtx )
                    v1 = cell(i)%vtx(k)
                    if (k == 1) then
                        v2 = cell(i)%vtx(2)
                        v3 = cell(i)%vtx(3)
                    elseif (k==2) then
                        v2 = cell(i)%vtx(1)
                        v3 = cell(i)%vtx(4)
                    elseif (k==3) then
                        v2 = cell(i)%vtx(2)
                        v3 = cell(i)%vtx(4)
                    elseif (k==4) then
                        v2 = cell(i)%vtx(1)
                        v3 = cell(i)%vtx(3)
                    end if
                    ! look for neighbor in the cells around v1
                    find_nghbr1 : do kk = 1, node(v1)%nc
                        ck = node(v1)%c(kk) ! identity of cell containing node kk
                        if (cell(ck)%nvtx == 4) then ! neighbor cell is a tet
                            vk1 = cell(ck)%vtx(1)
                            vk2 = cell(ck)%vtx(2)
                            vk3 = cell(ck)%vtx(3)
                            vk4 = cell(ck)%vtx(4)
                            if ( faceShareWithTet(v1,v2,v3,vk1,vk2,vk3,vk4) ) then
                                cell(i)%nnghbrs = cell(i)%nnghbrs + 1
                                temp(cell(i)%nnghbrs)    = ck ! Get the number of the neighbor cell
                                tempv(cell(i)%nnghbrs,1) = v1 ! Vertex 1  of the neighbor cell
                                tempv(cell(i)%nnghbrs,2) = v2 ! Vertex 1  of the neighbor cell
                                tempv(cell(i)%nnghbrs,3) = v3 ! Vertex 1  of the neighbor cell
                                tempc(cell(i)%nnghbrs)   = 3  ! number of common vertices w/ nghbr cell
                                exit find_nghbr1
                            end if
                        else if (cell(ck)%nvtx == 6) then
                            vk1 = cell(ck)%vtx(1)
                            vk2 = cell(ck)%vtx(2)
                            vk3 = cell(ck)%vtx(3)
                            vk4 = cell(ck)%vtx(4)
                            vk5 = cell(ck)%vtx(5)
                            vk6 = cell(ck)%vtx(6)                            
                            if ( faceShareWithPrism(v1,v2,v3,dummyScalar,vk1,vk2,vk3,vk4,vk5,vk6,'tri') ) then
                                cell(i)%nnghbrs = cell(i)%nnghbrs + 1
                                temp(cell(i)%nnghbrs)    = ck ! Get the number of the neighbor cell
                                tempv(cell(i)%nnghbrs,1) = v1 ! Vertex 1  of the neighbor cell
                                tempv(cell(i)%nnghbrs,2) = v2 ! Vertex 1  of the neighbor cell
                                tempv(cell(i)%nnghbrs,3) = v3 ! Vertex 1  of the neighbor cell
                                tempc(cell(i)%nnghbrs)   = 3  ! number of common vertices w/ nghbr cell
                                exit find_nghbr1
                            end if
                        else if (cell(ck)%nvtx == 8) then
                            cycle
                        else
                            write(*,*) "Cell ", i ," is not a tet, prism, or hex... Something is wrong."
                            stop
                        end if
                    end do find_nghbr1
                end do face_k_tet_loop
            else if (pnvtx == 6) then
                face_k_prism_loop : do k = 1, 5
                    if (k == 1) then
                        v1 = cell(i)%vtx(korder(k))
                        v2 = cell(i)%vtx(2)
                        v3 = cell(i)%vtx(3)
                    elseif ( k == 2 ) then
                        v3 = cell(i)%vtx(korder(k)) ! We need to revers the order of the top nodes
                        v2 = cell(i)%vtx(5)         ! That way we move CCW around the face (from 
                        v1 = cell(i)%vtx(6)         ! inside the cell).
                    elseif ( k == 3 ) then
                        v1 = cell(i)%vtx(korder(k))
                        v2 = cell(i)%vtx(3)
                        v3 = cell(i)%vtx(6)
                        v4 = cell(i)%vtx(4)
                    elseif ( k == 4 ) then
                        v1 = cell(i)%vtx(korder(k))
                        v2 = cell(i)%vtx(1)
                        v3 = cell(i)%vtx(4)
                        v4 = cell(i)%vtx(5)
                    elseif ( k == 5 ) then
                        v1 = cell(i)%vtx(korder(k))
                        v2 = cell(i)%vtx(2)
                        v3 = cell(i)%vtx(5)
                        v4 = cell(i)%vtx(6)
                    end if
                    if ((k == 1) .or. (k == 2)) then ! check the two tri faces first
                        ! Look for a neighbor in he cells around v1
                        find_nghbr2 : do kk = 1,node(v1)%nc ! for 1:number of cells connected to v1
                            ck = node(v1)%c(kk) ! Identity of cell containing node kk
                            if (cell(ck)%nvtx == 4) then ! neighbor cell is a tet
                                vk1 = cell(ck)%vtx(1)
                                vk2 = cell(ck)%vtx(2)
                                vk3 = cell(ck)%vtx(3)
                                vk4 = cell(ck)%vtx(4)
                                if ( faceShareWithTet(v1,v2,v3,vk1,vk2,vk3,vk4) ) then
                                    cell(i)%nnghbrs = cell(i)%nnghbrs + 1
                                    temp(cell(i)%nnghbrs) = ck ! Get the number of the neighbor cell
                                    tempv(cell(i)%nnghbrs,1) = v1 ! Vertex 1 of the neighbor cell
                                    tempv(cell(i)%nnghbrs,2) = v2 ! Vertex 2 of the neighbor cell
                                    tempv(cell(i)%nnghbrs,3) = v3 ! Vertex 3 of the neighbor cell
                                    tempc(cell(i)%nnghbrs)   = 3  ! number of common vertices w/ nghbr cell
                                    exit find_nghbr2
                                end if
                            elseif ( cell(ck)%nvtx == 6 ) then
                                vk1 = cell(ck)%vtx(1)
                                vk2 = cell(ck)%vtx(2)
                                vk3 = cell(ck)%vtx(3)
                                vk4 = cell(ck)%vtx(4)
                                vk5 = cell(ck)%vtx(5)
                                vk6 = cell(ck)%vtx(6)
                                if (faceShareWithPrism(v1, v2, v3, dummyScalar, vk1, vk2, vk3, vk4, vk5, vk6,'tri') ) then
                                    cell(i)%nnghbrs = cell(i)%nnghbrs + 1
                                    temp(cell(i)%nnghbrs) = ck ! Get the number of the neighbor cell
                                    tempv(cell(i)%nnghbrs,1) = v1 ! Vertex 1 of the neighbor cell
                                    tempv(cell(i)%nnghbrs,2) = v2 ! Vertex 2 of the neighbor cell
                                    tempv(cell(i)%nnghbrs,3) = v3 ! Vertex 3 of the neighbor cell
                                    tempc(cell(i)%nnghbrs)   = 3  ! number of common vertices w/ nghbr cell
                                    exit find_nghbr2
                                end if
                            else if ( cell(ck)%nvtx == 8 ) then
                                cycle
                            else 
                                write(*,*) "Cell ", i ," is not a tet, prism, or hex... Something is wrong."
                                stop                                
                            end if
                        end do find_nghbr2
                    else if (k == 3 .or. k == 4 .or. k ==5) then
                        find_nghbr3 : do kk = 1,node(v1)%nc ! for 1:number of cells connected to v1
                            ck = node(v1)%c(kk)
                            ! If the neighbor cell is a tet we can skip since we're looking at quads
                            if (cell(ck)%nvtx == 4) then
                                cycle
                            else if (cell(ck)%nvtx == 6) then ! neighbor is a prism
                                vk1 = cell(ck)%vtx(1)
                                vk2 = cell(ck)%vtx(2)
                                vk3 = cell(ck)%vtx(3)
                                vk4 = cell(ck)%vtx(4)
                                vk5 = cell(ck)%vtx(5)
                                vk6 = cell(ck)%vtx(6)
                                if ( faceShareWithPrism(v1,v2,v3,v4,vk1,vk2,vk3,vk4,vk5,vk6,'quad') ) then
                                    cell(i)%nnghbrs = cell(i)%nnghbrs + 1
                                    temp(cell(i)%nnghbrs) = ck ! Get the number of the neighbor cell
                                    tempv(cell(i)%nnghbrs,1) = v1 ! Vertex 1 of the neighbor cell
                                    tempv(cell(i)%nnghbrs,2) = v2 ! Vertex 2 of the neighbor cell
                                    tempv(cell(i)%nnghbrs,3) = v3 ! Vertex 3 of the neighbor cell
                                    tempv(cell(i)%nnghbrs,4) = v4 ! Vertex 4 of the neighbor cell
                                    tempc(cell(i)%nnghbrs)   = 4  ! number of common vertices w/ nghbr cell
                                    exit find_nghbr3
                                end if  
                            else if ( cell(ck)%nvtx == 8 ) then
                                vk1 = cell(ck)%vtx(1)
                                vk2 = cell(ck)%vtx(2)
                                vk3 = cell(ck)%vtx(3)
                                vk4 = cell(ck)%vtx(4)
                                vk5 = cell(ck)%vtx(5)
                                vk6 = cell(ck)%vtx(6)
                                vk7 = cell(ck)%vtx(7)
                                vk8 = cell(ck)%vtx(8)
                                if ( faceShareWithHex(v1,v2,v3,v4,vk1,vk2,vk3,vk4,vk5,vk6,vk7,vk8) ) then
                                    cell(i)%nnghbrs = cell(i)%nnghbrs + 1
                                    temp(cell(i)%nnghbrs) = ck ! Get the number of the neighbor cell
                                    tempv(cell(i)%nnghbrs,1) = v1 ! Vertex 1 of the neighbor cell
                                    tempv(cell(i)%nnghbrs,2) = v2 ! Vertex 2 of the neighbor cell
                                    tempv(cell(i)%nnghbrs,3) = v3 ! Vertex 3 of the neighbor cell
                                    tempv(cell(i)%nnghbrs,4) = v4 ! Vertex 4 of the neighbor cell
                                    tempc(cell(i)%nnghbrs)   = 4  ! number of common vertices w/ nghbr cell
                                    exit find_nghbr3
                                end if  
                            else
                                write(*,*) "Cell ", i ," is not a tet, prism, or hex... Something is wrong."
                                stop                
                            end if
                        end do find_nghbr3
                    end if
                end do face_k_prism_loop
            else if (pnvtx == 8) then
                face_k_hex_loop : do k = 1, 6
                    if ( k == 1) then
                        v1 = cell(i)%vtx(1)
                        v2 = cell(i)%vtx(2)
                        v3 = cell(i)%vtx(3)
                        v4 = cell(i)%vtx(4)
                    elseif ( k == 2) then
                        v1 = cell(i)%vtx(5)
                        v2 = cell(i)%vtx(6)
                        v3 = cell(i)%vtx(2)
                        v4 = cell(i)%vtx(1)
                    elseif ( k == 3) then
                        v1 = cell(i)%vtx(6)
                        v2 = cell(i)%vtx(7)
                        v3 = cell(i)%vtx(3)
                        v4 = cell(i)%vtx(2)
                    elseif ( k == 4) then
                        v1 = cell(i)%vtx(7)
                        v2 = cell(i)%vtx(8)
                        v3 = cell(i)%vtx(4)
                        v4 = cell(i)%vtx(3)
                    elseif ( k == 5) then
                        v1 = cell(i)%vtx(8)
                        v2 = cell(i)%vtx(5)
                        v3 = cell(i)%vtx(1)
                        v4 = cell(i)%vtx(4)
                    elseif ( k == 6) then
                        v1 = cell(i)%vtx(8)
                        v2 = cell(i)%vtx(7)
                        v3 = cell(i)%vtx(6)
                        v4 = cell(i)%vtx(5)
                    endif
                    find_nghbr4 : do kk = 1,node(v1)%nc ! for 1:number of cells connected to v1
                            ck = node(v1)%c(kk)
                            ! If the neighbor cell is a tet we can skip since we're looking at quads
                            if (cell(ck)%nvtx == 4) then
                                cycle
                            else if (cell(ck)%nvtx == 6) then ! neighbor is a prism
                                vk1 = cell(ck)%vtx(1)
                                vk2 = cell(ck)%vtx(2)
                                vk3 = cell(ck)%vtx(3)
                                vk4 = cell(ck)%vtx(4)
                                vk5 = cell(ck)%vtx(5)
                                vk6 = cell(ck)%vtx(6)
                                if ( faceShareWithPrism(v1,v2,v3,v4,vk1,vk2,vk3,vk4,vk5,vk6,'quad') ) then
                                    cell(i)%nnghbrs = cell(i)%nnghbrs + 1
                                    temp(cell(i)%nnghbrs) = ck ! Get the number of the neighbor cell
                                    tempv(cell(i)%nnghbrs,1) = v1 ! Vertex 1 of the neighbor cell
                                    tempv(cell(i)%nnghbrs,2) = v2 ! Vertex 2 of the neighbor cell
                                    tempv(cell(i)%nnghbrs,3) = v3 ! Vertex 3 of the neighbor cell
                                    tempv(cell(i)%nnghbrs,4) = v4 ! Vertex 4 of the neighbor cell
                                    tempc(cell(i)%nnghbrs)   = 4  ! number of common vertices w/ nghbr cell
                                    exit find_nghbr4
                                end if  
                            else if ( cell(ck)%nvtx == 8 ) then
                                vk1 = cell(ck)%vtx(1)
                                vk2 = cell(ck)%vtx(2)
                                vk3 = cell(ck)%vtx(3)
                                vk4 = cell(ck)%vtx(4)
                                vk5 = cell(ck)%vtx(5)
                                vk6 = cell(ck)%vtx(6)
                                vk7 = cell(ck)%vtx(7)
                                vk8 = cell(ck)%vtx(8)
                                if ( faceShareWithHex(v1,v2,v3,v4,vk1,vk2,vk3,vk4,vk5,vk6,vk7,vk8) ) then
                                    cell(i)%nnghbrs = cell(i)%nnghbrs + 1
                                    temp(cell(i)%nnghbrs) = ck ! Get the number of the neighbor cell
                                    tempv(cell(i)%nnghbrs,1) = v1 ! Vertex 1 of the neighbor cell
                                    tempv(cell(i)%nnghbrs,2) = v2 ! Vertex 2 of the neighbor cell
                                    tempv(cell(i)%nnghbrs,3) = v3 ! Vertex 3 of the neighbor cell
                                    tempv(cell(i)%nnghbrs,4) = v4 ! Vertex 4 of the neighbor cell
                                    tempc(cell(i)%nnghbrs)   = 4  ! number of common vertices w/ nghbr cell
                                    exit find_nghbr4
                                end if  
                            else
                                write(*,*) "Cell ", i ," is not a tet, prism, or hex... Something is wrong."
                                stop                
                            end if
                        end do find_nghbr4
                        
                        
                end do face_k_hex_loop
            end if cell_type
            ! Store the nighbors in the neghbor list
            allocate( cell(i)%nghbr(   cell(i)%nnghbrs )   )
            allocate( cell(i)%ncommonv(cell(i)%nnghbrs )   )
            allocate( cell(i)%commonv( cell(i)%nnghbrs,4 ) )
            ! Allocate the jacobian off diagonal term while were at it
            allocate(  jac(i)%off_diag(5,5,cell(i)%nnghbrs))
            do k = 1,cell(i)%nnghbrs
                cell(i)%nghbr(k) = temp(k)              ! kth neighbor cell
                cell(i)%ncommonv(k) = tempc(k)
                do j = 1,cell(i)%ncommonv(k)
                    cell(i)%commonv(k,j) = tempv(k,j)
                end do
            end do
        end do cell_assoc_loop

        !---------------------------------------------------------------
        !---------------------------------------------------------------
        ! (4) Construction of face data.
        !---------------------------------------------------------------
        
        ! Note: The face data contains only the interior faces, which are shared by two cells
        ! Note: Define the face as pointing from cell1 to cell2 where cell2 > cell1
        !       (i.e. smaller cell number to larger cell number).
        !       So, a face is oriented from a cell number to a larger cell number.
        write(*,*)
        write(*,*) " Constructing face data"
        write(*,*)

        !---------------------------------------
        ! First count the number of total interior faces
        nfaces = 0
        do i = 1,ncells
            do k = 1,cell(i)%nnghbrs
                ! This identifies a face and avoids a double count.
                ! Ex. if cell 5 and cell 6 share a face, when we get to i = 5 that 
                ! face will be added to the count.  However, when i=6 and we loop 
                ! back over the same face it won't be counted since 6<5 = .false.
                if (i < cell(i)%nghbr(k)) then
                    nfaces = nfaces + 1
                end if
            end do
        end do
        !---------------------------------------
        ! Allocate and fill the face arrays.
        ! Switched to column major (Matlab code is row major, which is gonna make copying stuff fun)
        allocate(face(         7,nfaces)) ! face list: (/Left_cell, Right_cell,#_of_sides, node1:node_end/)
        allocate(face_nrml(    3,nfaces)) ! unit face normal vector
        allocate(face_nrml_mag(  nfaces)) ! magnitude of face vector
        allocate(face_centroid(3,nfaces)) ! location of face_centroid
        ! We need to re-count nfaces, so, initialize it again
        nfaces = 0
        ! Now, construct the face arrays
        cell_loop_f : do i = 1,ncells
            neighbor_loop : do k = 1,cell(i)%nnghbrs
                if (i<cell(i)%nghbr(k)) then
                    ! Face is found when the neighbor cell nunmber is greater than i
                    nfaces = nfaces + 1
                    !           Left(2)
                    !        o---o---------o
                    !       .    .          .
                    !      .     .           .
                    !     .      .normal      .
                    !    .  Left .--->  Right  .
                    !   .    i   .       k      .
                    !  .         .               .
                    ! o----------o----------------o
                    !          Right(1)
                    !
                    ! Left Cell
                    face(1,nfaces) = i
                    ! Right Cell
                    face(2,nfaces) = cell(i)%nghbr(k)
                    ! Number of points
                    face(3,nfaces) = cell(i)%ncommonv(k)
                    do j = 1,face(3, nfaces)
                        face(j+3, nfaces) = cell(i)%commonv(k,j)
                    end do
                    ! Root vertex (seen from i)
                    x1 = x(cell(i)%commonv(k,1))
                    y1 = y(cell(i)%commonv(k,1))
                    z1 = z(cell(i)%commonv(k,1))
                    ! Root vertex (seen from i)
                    x2 = x(cell(i)%commonv(k,2))
                    y2 = y(cell(i)%commonv(k,2))
                    z2 = z(cell(i)%commonv(k,2))
                    ! Root vertex (seen from i)
                    x3 = x(cell(i)%commonv(k,3))
                    y3 = y(cell(i)%commonv(k,3))
                    z3 = z(cell(i)%commonv(k,3))
                    if (face(3,nfaces) == 4) then ! quad face
                        ! Root vertex (seen from i)
                        x4 = x(cell(i)%commonv(k,4))
                        y4 = y(cell(i)%commonv(k,4))
                        z4 = z(cell(i)%commonv(k,4))
                        call quadNormal(x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4,face_nrml(1:3,nfaces),face_nrml_mag(nfaces))
                        face_nrml(1:3,nfaces) = -one*face_nrml(1:3,nfaces) ! ensure the vector is pointing out
                        call get_quad_face_centroid(x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4, face_centroid(1:3,nfaces))
                    else if (face(3,nfaces) == 3) then
                        call triNormal(x1,y1,z1, x2,y2,z2, x3,y3,z3, face_nrml(1:3,nfaces),face_nrml_mag(nfaces))
                        face_nrml(1:3,nfaces) = -one*face_nrml(1:3,nfaces) ! ensure the vector is pointing out
                        call get_tri_face_centroid(x1,y1,z1, x2,y2,z2, x3,y3,z3, face_centroid(1:3,nfaces))
                    else
                        write(*,*) " Wrong face size. STOP!"
                        stop
                    end if
                end if
            end do neighbor_loop
        end do cell_loop_f

        ! Create kth_nghbr arrays
        allocate(kth_nghbr_of_1(nfaces))
        allocate(kth_nghbr_of_2(nfaces))
        
        face_nghbr_loop : do i = 1,nfaces
            c1 = face(1,i)
            c2 = face(2,i)
            ! loop over c1 neighbors to find c2
            do k = 1,cell(c1)%nnghbrs
                if ( c2 == cell(c1)%nghbr(k)) then
                    kth_nghbr_of_1(i) = k ! c2 is the kth neighbor of c1
                end if
            end do
            ! repeat for cell 2
            do k = 1,cell(c2)%nnghbrs
                if ( c1 == cell(c2)%nghbr(k)) then
                    kth_nghbr_of_2(i) = k
                end if
            end do
        end do face_nghbr_loop


        !-------------------------------------------------------------------------------
        ! Compute and store the boundary face normals.
        !-------------------------------------------------------------------------------

        ! Boundary face j consists of nodes j and j+1.
        !
        !  Interior domain      /
        !                      /
        !              /\     o
        !             /  \   /
        !            / ck \ /   Outside the domain
        ! --o-------o------o
        !           j   |  j+1
        !               |   
        !               v Face normal for the face j.
        !
        ! ck = bcell, the cell having the boundary face j.
        !


        ! Note: in this case the boundary faces are all of the tria and quad 
        write(*,*)
        write(*,*) " Compute boundary face normals...."
        write(*,*)
        
        allocate(bound(nb))
        ! count and sort boundary face in each zone
        ! Initialize nbfaces
        boundary_parts : do i = 1,nb
            bound(i)%nbfaces = 0
        end do boundary_parts
        ! Count bfaces for each zone
        do i = 1,ntria ! first sort tria
            bound(tria(i,4))%nbfaces = bound(tria(i,4))%nbfaces + 1
        end do
        do i = 1,nquad ! next sort quads
            bound(quad(i,5))%nbfaces = bound(quad(i,5))%nbfaces + 1
        end do
        ! allocate bfaces and re-initialize nbfaces
        do i = 1,nb
            allocate(bound(i)%bfaces(5,bound(i)%nbfaces)) ! (/#_of_nodes v1:v_end/)
            bound(i)%nbfaces = 0
        end do
        ! Sort bfaces
        do i = 1,ntria
            bound(tria(i,4))%nbfaces = bound(tria(i,4))%nbfaces + 1
            bound(tria(i,4))%bfaces(1:4,bound(tria(i,4))%nbfaces) = (/3, tria(i,1), tria(i,2), tria(i,3)/)
        end do
        do i = 1,nquad ! next sort quads
            bound(quad(i,5))%nbfaces = bound(quad(i,5))%nbfaces + 1
            bound(quad(i,5))%bfaces(1:5,bound(quad(i,5))%nbfaces) = (/4, quad(i,1), quad(i,2), quad(i,3), quad(i,4)/)
        end do
        ! Loop through each face to get its normal
        ! TODO: implement bface_centroid to this loop (done. i think...)
        do i = 1,nb
            allocate(bound(i)%bface_nrml(    3,bound(i)%nbfaces))
            allocate(bound(i)%bface_center(  3,bound(i)%nbfaces))
            allocate(bound(i)%bface_nrml_mag(bound(i)%nbfaces  ))
            do j = 1,bound(i)%nbfaces
                if ( bound(i)%bfaces(1,j) ==3 ) then ! tri face
                    v1 = bound(i)%bfaces(4,j) ! flipped for consistency
                    v2 = bound(i)%bfaces(3,j)
                    v3 = bound(i)%bfaces(2,j)
                else if ( bound(i)%bfaces(1,j) == 4 ) then ! quad face
                    v1 = bound(i)%bfaces(5,j) ! flipped for consistency
                    v2 = bound(i)%bfaces(4,j)
                    v3 = bound(i)%bfaces(3,j)
                    v4 = bound(i)%bfaces(2,j)
                end if
                ! Root vertex (seen from i)
                x1 = x(v1)
                y1 = y(v1)
                z1 = z(v1)
                ! Root vertex (seen from i)
                x2 = x(v2)
                y2 = y(v2)
                z2 = z(v2)
                ! Root vertex (seen from i)
                x3 = x(v3)
                y3 = y(v3)
                z3 = z(v3)
                if (bound(i)%bfaces(1,j) == 4) then ! quad face
                    ! Root vertex (seen from i)
                    x4 = x(v4)
                    y4 = y(v4)
                    z4 = z(v4)
                    call quadNormal(x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4, &
                                    bound(i)%bface_nrml(1:3,j), bound(i)%bface_nrml_mag(j))
                    call get_quad_face_centroid(x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4, &
                            bound(i)%bface_center(1:3,j))
                else if (bound(i)%bfaces(1,j) == 3) then ! tri face
                    call  triNormal(x1,y1,z1, x2,y2,z2, x3,y3,z3, &
                                    bound(i)%bface_nrml(1:3,j), bound(i)%bface_nrml_mag(j))
                    ! Note doesn't need to be flipped since we're looking outside in
                    call get_tri_face_centroid(x1,y1,z1, x2,y2,z2, x3,y3,z3, &
                                    bound(i)%bface_center(1:3,j))
                end if
            end do
        end do
        
        ! Loop through each face to get its attached cell and face normal
        do i = 1,nb
            allocate(bound(i)%bcell(bound(i)%nbfaces))
            do j = 1,bound(i)%nbfaces
                ! get the vertices again
                if ( bound(i)%bfaces(1,j) ==3 ) then ! tri face
                    v1 = bound(i)%bfaces(4,j) ! flipped for consistency
                    v2 = bound(i)%bfaces(3,j)
                    v3 = bound(i)%bfaces(2,j)
                    find_bnghbr : do kk = 1,node(v1)%nc
                        ck = node(v1)%c(kk)
                        if (cell(ck)%nvtx == 4) then
                            vk1 = cell(ck)%vtx(1)
                            vk2 = cell(ck)%vtx(2)
                            vk3 = cell(ck)%vtx(3)
                            vk4 = cell(ck)%vtx(4)
                            if (faceShareWithTet(v1,v2,v3, vk1,vk2,vk3,vk4)) then
                                bound(i)%bcell(j) = ck
                                exit find_bnghbr
                            end if
                        else if ( cell(ck)%nvtx == 6 ) then ! prism
                            vk1 = cell(ck)%vtx(1)
                            vk2 = cell(ck)%vtx(2)
                            vk3 = cell(ck)%vtx(3)
                            vk4 = cell(ck)%vtx(4)
                            vk5 = cell(ck)%vtx(5)
                            vk6 = cell(ck)%vtx(6)
                            if ( faceShareWithPrism(v1,v2,v3,dummyScalar, &
                                        vk1,vk2,vk3,vk4,vk5,vk6,'tri')) then
                                bound(i)%bcell(j) = ck
                                exit find_bnghbr
                            end if   
                        elseif ( cell(ck)%nvtx == 8 ) then ! hex
                            cycle find_bnghbr
                        end if
                    end do find_bnghbr
                else if (bound(i)%bfaces(1,j) == 4) then ! quad face
                    v1 = bound(i)%bfaces(5,j) ! flipped for consistency
                    v2 = bound(i)%bfaces(4,j)
                    v3 = bound(i)%bfaces(3,j)
                    v4 = bound(i)%bfaces(2,j)
                    find_qnghbr : do kk = 1,node(v1)%nc
                        ck = node(v1)%c(kk)
                        if ( cell(ck)%nvtx == 4 ) then
                            cycle find_qnghbr
                        else if ( cell(ck)%nvtx == 6 ) then ! prism
                            vk1 = cell(ck)%vtx(1)
                            vk2 = cell(ck)%vtx(2)
                            vk3 = cell(ck)%vtx(3)
                            vk4 = cell(ck)%vtx(4)
                            vk5 = cell(ck)%vtx(5)
                            vk6 = cell(ck)%vtx(6)
                            if ( faceShareWithPrism(v1,v2,v3,v4, &
                                        vk1,vk2,vk3,vk4,vk5,vk6,'quad')) then
                                bound(i)%bcell(j) = ck
                                exit find_qnghbr
                            end if
                        elseif ( cell(ck)%nvtx == 8 ) then ! hex
                            vk1 = cell(ck)%vtx(1)
                            vk2 = cell(ck)%vtx(2)
                            vk3 = cell(ck)%vtx(3)
                            vk4 = cell(ck)%vtx(4)
                            vk5 = cell(ck)%vtx(5)
                            vk6 = cell(ck)%vtx(6)
                            vk7 = cell(ck)%vtx(7)
                            vk8 = cell(ck)%vtx(8)
                            if ( faceShareWithHex(v1,v2,v3,v4,vk1,vk2,vk3,vk4,vk5,vk6,vk7,vk8) ) then
                                bound(i)%bcell(j) = ck
                                exit find_qnghbr
                            end if  
                        end if
                    end do find_qnghbr
                end if
            end do
        end do

        ! Now it's time to compute the cell volume
        ! this will be done using the divergence theorem.  Method has been
        ! developed from Fluid Mechanics 101 youtube video entitled:
        ! "Calculating the Cell Volume"
        ! Vol_cell = sum((1/3)(x_f*n_f)A_f)
        ! x_f = vector from the cell center to the face center 
        ! n_f = unit face normal
        ! A_f = face area
        !              o
        !             / \
        !            /   \
        !           /     \
        !          /       \
        !         /         \
        !        /	         \
        !       |	          |
        !   n_f |	          |
        ! <-----|<-----o      |
        !       |  x_f        |
        !       |_____________|

        ! Initialize cell volumes to zero
        do i = 1,ncells
            cell(i)%vol = 0
        end do
        ! first loop through internal faces
        do i = 1,nfaces
            c1 = face(1,i)
            c2 = face(2,i)
            x1 = cell(c1)%xc
            y1 = cell(c1)%yc
            z1 = cell(c1)%zc
            x2 = cell(c2)%xc
            y2 = cell(c2)%yc
            z2 = cell(c2)%zc
            x_f = face_centroid(1:3,i) - (/x1, y1, z1/)
            dot_result1 = dot_product(x_f,  face_nrml(1:3,i)) ! face normal is out of cell 1
            x_f = face_centroid(1:3,i) - (/x2, y2, z2/)
            dot_result2 = dot_product(x_f, -face_nrml(1:3,i)) ! face normal is into cell 2
            cell(c1)%vol = (one/three) * dot_result1 * face_nrml_mag(i) + cell(c1)%vol
            cell(c2)%vol = (one/three) * dot_result2 * face_nrml_mag(i) + cell(c2)%vol
        end do

        ! Now boundary faces
        do i = 1,nb
            do j = 1,bound(i)%nbfaces
                c1 = bound(i)%bcell(j)
                x1 = cell(c1)%xc
                y1 = cell(c1)%yc
                z1 = cell(c1)%zc
                x_f = bound(i)%bface_center(1:3,j) - (/x1, y1, z1/)
                dot_result1 = dot_product(x_f,bound(i)%bface_nrml(1:3,j))
                cell(c1)%vol = (one/three) * dot_result1 * bound(i)%bface_nrml_mag(j) + cell(c1)%vol
            end do
        end do

        ! NOTE: I still need to deallocate the node array!!
        deallocate(node) ! I don't wanna play with you anymore...
        ! nevermind...

        !allocate( local_cell(ncells))
        !local_cell = cell
        !allocate(local_face(7,nfaces))
        !local_face  = face
        !allocate(local_face_nrml(    3,nfaces)) ! unit face normal vector
        !allocate(local_face_nrml_mag(  nfaces)) ! magnitude of face vector
        !allocate(local_face_centroid(3,nfaces)) ! location of face_centroid
        !local_face_nrml = face_nrml
        !local_face_nrml_mag = face_nrml_mag
        !local_face_centroid = face_centroid
        !allocate(local_bound(nb))
        !local_bound = bound
        !
        ! Now to do some mesh checks
        write(*,*)
        write(*,*)
        write(*,*) "  --- Checking cell centroids --- "
        do i = 1,ncells
            if (cell(i)%nvtx == 4) then
                v1 = cell(i)%vtx(1)
                v2 = cell(i)%vtx(2)
                v3 = cell(i)%vtx(3)
                v4 = cell(i)%vtx(4)
                inside = checkInsideTet(x(v1),y(v1),z(v1), x(v2),y(v2),z(v2), &
                        x(v3),y(v3),z(v3), x(v4),y(v4),z(v4), &
                        cell(i)%xc,cell(i)%yc,cell(i)%zc)
            else if (cell(i)%nvtx == 6) then
                inside = .true. ! not implemented yet...
            else if (cell(i)%nvtx == 8) then
                inside = .true. ! not implemented yet...
            end if
            if (inside .eqv. .false.) then
                write(*,*) " Cell ", i ," centroid is outside of cell.  Cannot continue. Stop!"
                stop
            end if
        end do
        write(*,*) " Check complete!"
        write(*,*)

        ! Skipping volume check for now, use the green-gauss theorem as above with the boundary faces

        ! Boundary Normal Check
        machine0_1 = epsilon(one)
        write(*,*) "     Machine zero relative to 1 = ",  machine0_1
        sum_machine0_1 = zero
        do i = 1, nb
            do j = 1, bound(i)%nbfaces
                sum_machine0_1  = sum_machine0_1  + bound(i)%bface_nrml_mag(j)*machine0_1 
            end do
        end do
        write(*,*) "     Machine zero relative to sum of boundary face mag = ", sum_machine0_1
        write(*,*)
        sum_bnormal = (/zero, zero, zero/)
        do i = 1,nb
            do j = 1,bound(i)%nbfaces
                sum_bnormal(1) = sum_bnormal(1) + bound(i)%bface_nrml(1,j) * bound(i)%bface_nrml_mag(j)
                sum_bnormal(2) = sum_bnormal(2) + bound(i)%bface_nrml(2,j) * bound(i)%bface_nrml_mag(j)
                sum_bnormal(3) = sum_bnormal(3) + bound(i)%bface_nrml(3,j) * bound(i)%bface_nrml_mag(j)
            end do
        end do
        write(*,*) "     Sum of boundary face normal (nx) = ", sum_bnormal(ix)
        write(*,*) "     Sum of boundary face normal (ny) = ", sum_bnormal(iy)
        write(*,*) "     Sum of boundary face normal (nz) = ", sum_bnormal(iz)

        if ( abs(sum_bnormal(ix)) > 50.0_p2*sum_machine0_1 .or. &
             abs(sum_bnormal(iy)) > 50.0_p2*sum_machine0_1 .or. &
             abs(sum_bnormal(iz)) > 50.0_p2*sum_machine0_1      ) then
            write(*,*) " Boundary face vector sum is larger than machine zero....."
            write(*,*) " Something is wrong. Stop."
            stop
        end if
        !--------------------------------------------------------
        ! Face normal sum check for each cell
        !--------------------------------------------------------
        write(*,*)
        write(*,*)
        write(*,*) " --- Check the face-normal sum over faces of a cell:"
        write(*,*)

        allocate(sum_face_normal(3,ncells))
        sum_face_normal = zero

        !----------------------------------------------------
        ! Accumulate face normals at cells by looping over interior faces.

        do i = 1, nfaces
            cell1 = face(1,i)
            sum_face_normal(ix,cell1) = sum_face_normal(ix,cell1) + face_nrml(ix,i)*face_nrml_mag(i)
            sum_face_normal(iy,cell1) = sum_face_normal(iy,cell1) + face_nrml(iy,i)*face_nrml_mag(i)
            sum_face_normal(iz,cell1) = sum_face_normal(iz,cell1) + face_nrml(iz,i)*face_nrml_mag(i)
            cell2 = face(2,i)
            sum_face_normal(ix,cell2) = sum_face_normal(ix,cell2) - face_nrml(ix,i)*face_nrml_mag(i)
            sum_face_normal(iy,cell2) = sum_face_normal(iy,cell2) - face_nrml(iy,i)*face_nrml_mag(i)
            sum_face_normal(iz,cell2) = sum_face_normal(iz,cell2) - face_nrml(iz,i)*face_nrml_mag(i)
        end do

        !----------------------------------------------------
        ! Add boundary face normal contributions to cells by looping over boundary faces.

        ! Loop over boundary segments
        do i = 1,nb
            do j = 1,bound(i)%nbfaces
                cell1 = bound(i)%bcell(j)
                sum_face_normal(ix,cell1) = sum_face_normal(ix,cell1) + bound(i)%bface_nrml(ix,j)*bound(i)%bface_nrml_mag(j)
                sum_face_normal(iy,cell1) = sum_face_normal(iy,cell1) + bound(i)%bface_nrml(iy,j)*bound(i)%bface_nrml_mag(j)
                sum_face_normal(iz,cell1) = sum_face_normal(iz,cell1) + bound(i)%bface_nrml(iz,j)*bound(i)%bface_nrml_mag(j)
            end do
        end do
        write(*,*) "   Display the maximum over all cells (must be zero):"
        write(*,*) "   (the sum of face normals over each cell)"
        write(*,*)
        write(*,*) "     Max of |sum_faces face normal (nx)| = ", maxval(abs(sum_face_normal(ix,:)))
        write(*,*) "     Max of |sum_faces face normal (ny)| = ", maxval(abs(sum_face_normal(iy,:)))
        write(*,*) "     Max of |sum_faces face normal (nz)| = ", maxval(abs(sum_face_normal(iz,:)))
        
        
        !----------------------------------------------------
        ! Check the maximum sum over all cells. Must be zero.

        !Use the maximum sqrt(vol) as a reference magnitude for checking zero face normal sum.

        ref_mag = sqrt( maxval( cell(1:ncells)%vol) )



        write(*,*)
        write(*,*) "              Reference magnitude = ", ref_mag
        write(*,*) "      Machine zero w.r.t. ref mag = ", ref_mag*machine0_1
        write(*,*)
        do i = 1,ncells
            if ( (abs(sum_face_normal( ix, i ))) > 50.0_p2*ref_mag*machine0_1 .or. &
             (abs(sum_face_normal( iy, i ))) > 50.0_p2*ref_mag*machine0_1 .or. &
             (abs(sum_face_normal( iz, i ))) > 50.0_p2*ref_mag*machine0_1       ) then
                write(*,*) "Sum face noraml_x(", i ,") = ",sum_face_normal(ix,i)
                write(*,*) "Sum face noraml_y(", i ,") = ",sum_face_normal(iy,i)
                write(*,*) "Sum face noraml_z(", i ,") = ",sum_face_normal(iz,i)
                write(*,*)
            end if
        end do
        if ( maxval(abs(sum_face_normal( ix, 1:ncells ))) > 50.0_p2*ref_mag*machine0_1 .or. &
             maxval(abs(sum_face_normal( iy, 1:ncells ))) > 50.0_p2*ref_mag*machine0_1 .or. &
             maxval(abs(sum_face_normal( iz, 1:ncells ))) > 50.0_p2*ref_mag*machine0_1       ) then

            write(*,*) " Max face vector sum over a cell is larger than machine zero... Something is wrong. Stop."

            stop

        endif

        deallocate(sum_face_normal)

        !Any other check?

        write(*,*)
        write(*,*)
        write(*,*) " Finished Verifying the ccfv grid data..."

        !------------------------------------------------------------------------------------
        ! Compute mesh spacing statistics.
        !
        !   heffn     = Effecrtive spacing based on # of nodes
        !   heffc     = Effecrtive spacing based on # of cells
        !   heffv     = Average of sqrt(volume).
        !   heffv_min = Minimum sqrt(volume).
        !   heffv_max = Maximum sqrt(volume).
        !

        write(*,*)
        write(*,*)
        write(*,*) " Compute effective mesh spacings"
        write(*,*)

        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        heffn     = sqrt( one/real(nnodes,p2) ) !Effecrtive spacing based on # of nodes
        heffc     = sqrt( one/real(ncells,p2) ) !Effecrtive spacing based on # of cells

        !Compute average, min, and max effective mesh spacing based on sqrt(volume).

        heffv     = sqrt( cell(1)%vol )
        heffv_min = sqrt( cell(1)%vol )
        heffv_max = sqrt( cell(1)%vol )

        do i = 2, ncells
            heffv     = heffv + sqrt( cell(i)%vol )
            heffv_min = min( heffv_min, sqrt( cell(i)%vol ) )
            heffv_max = max( heffv_max, sqrt( cell(i)%vol ) )
        end do

        heffv     = heffv / real(ncells,p2)

        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        write(*,*) "  heffn     = ", heffn
        write(*,*) "  heffc     = ", heffc
        write(*,*)
        write(*,*) "  heffv     = ", heffv
        write(*,*) "  heffv_min = ", heffv_min
        write(*,*) "  heffv_max = ", heffv_max
        write(*,*)
        write(*,*)

        !------------------------------------------------------------------------------------
        !------------------------------------------------------------------------------------
        !------------------------------------------------------------------------------------

        write(*,*)
        write(*,*) " End of Construct ccfv grid data and verifty them." 
        write(*,*) "-------------------------------------------------------"
        write(*,*) "-------------------------------------------------------"
        write(*,*) "-------------------------------------------------------"
        write(*,*)

        ! We're done!  I'm tired.  I hope it works!
    end subroutine construct_ccfv_grid_data

    logical function faceShareWithTet(f1,f2,f3,c1,c2,c3,c4)
        implicit none
        integer, intent(in) ::  f1,  f2,  f3       ! Vertices of host cell face
        integer, intent(in) ::  c1,  c2,  c3,  c4  ! Vertices of candidate neighbor cells
        faceShareWithTet = .false.
        ! face = [f1, f2, f3]
        ! cell = [c1, c2, c3, c4]

        ! a in bottom right corner
        if ((f1 == c2 .and. f2 == c1 .and. f3 == c3) .or. &
            (f1 == c3 .and. f2 == c1 .and. f3 == c4) .or. &
            (f1 == c4 .and. f2 == c1 .and. f3 == c2)) then
                faceShareWithTet = .true.
                return
        end if
        ! b in bottom right corner
        if ((f1 == c3 .and. f2 == c2 .and. f3 == c1) .or. &
            (f1 == c1 .and. f2 == c2 .and. f3 == c4) .or. &
            (f1 == c4 .and. f2 == c2 .and. f3 == c3)) then
                faceShareWithTet = .true.
                return
        end if
        ! c in bottom right corner
        if ((f1 == c2 .and. f2 == c3 .and. f3 == c4) .or. &
            (f1 == c1 .and. f2 == c3 .and. f3 == c2) .or. &
            (f1 == c4 .and. f2 == c3 .and. f3 == c1)) then
                faceShareWithTet = .true.
                return
        end if
        ! d in bottom right corner
        if ((f1 == c1 .and. f2 == c4 .and. f3 == c3) .or. &
            (f1 == c2 .and. f2 == c4 .and. f3 == c1) .or. &
            (f1 == c3 .and. f2 == c4 .and. f3 == c2)) then
                faceShareWithTet = .true.
                return
        end if        
    end function faceShareWithTet
    logical function faceShareWithPrism(f1,f2,f3,f4,c1,c2,c3,c4,c5,c6,faceType)
        implicit none
        integer, intent(in) ::  f1,  f2,  f3,  f4            ! Vertices of host cell face
        integer, intent(in) ::  c1,  c2,  c3,  c4,  c5,  c6  ! Vertices of candidate neighbor cells
        character(len=*), intent(in) :: faceType
        faceShareWithPrism = .false.
        ! check settings are correct
        if (.not.((faceType == 'tri' .and. f4 == 0) .or. (faceType == 'quad'))) then
            write(*,*) " wrong setting prism assoc."
            stop
        end if
        if (faceType == 'tri') then
            ! test if match with bottom (c1 c2 c3)
            if ((f1 == c1 .and. f2 == c3 .and. f3 == c2) .or. &
                (f1 == c3 .and. f2 == c2 .and. f3 == c1) .or. &
                (f1 == c2 .and. f2 == c1 .and. f3 == c3)) then
                    faceShareWithPrism = .true.
                    return
            end if
            ! Test if match with top face (vv4 v5 v6)
            if ((f1 == c6 .and. f2 == c4 .and. f3 == c5) .or. &
                (f1 == c4 .and. f2 == c5 .and. f3 == c6) .or. &
                (f1 == c5 .and. f2 == c6 .and. f3 == c4)) then
                    faceShareWithPrism = .true.
                    return
            end if
        elseif (faceType == 'quad') then
            ! test if match with side 1
            if ((f1 == c3 .and. f2 == c1 .and. f3 == c4 .and. f4 == c6) .or. &
                (f1 == c1 .and. f2 == c4 .and. f3 == c6 .and. f4 == c3) .or. &
                (f1 == c4 .and. f2 == c6 .and. f3 == c3 .and. f4 == c1) .or. &
                (f1 == c6 .and. f2 == c3 .and. f3 == c1 .and. f4 == c4)) then
                faceShareWithPrism = .true.
                return
            end if
            ! test if match with side 2
            if ((f1 == c1 .and. f2 == c2 .and. f3 == c5 .and. f4 == c4) .or. &
                (f1 == c2 .and. f2 == c5 .and. f3 == c4 .and. f4 == c1) .or. &
                (f1 == c5 .and. f2 == c4 .and. f3 == c1 .and. f4 == c2) .or. &
                (f1 == c4 .and. f2 == c1 .and. f3 == c2 .and. f4 == c5)) then
                faceShareWithPrism = .true.
                return
            end if
            ! test if match with side 3
            if ((f1 == c2 .and. f2 == c3 .and. f3 == c6 .and. f4 == c5) .or. &
                (f1 == c3 .and. f2 == c6 .and. f3 == c5 .and. f4 == c2) .or. &
                (f1 == c6 .and. f2 == c5 .and. f3 == c2 .and. f4 == c3) .or. &
                (f1 == c5 .and. f2 == c2 .and. f3 == c3 .and. f4 == c6)) then
                faceShareWithPrism = .true.
                return
            end if
        end if
    end function faceShareWithPrism

    logical function faceShareWithhex(f1,f2,f3,f4,c1,c2,c3,c4,c5,c6,c7,c8)
        implicit none
        integer, intent(in) ::  f1,  f2,  f3,  f4                    ! Vertices of host cell face
        integer, intent(in) ::  c1,  c2,  c3,  c4,  c5,  c6, c7, c8  ! Vertices of candidate neighbor cells
        
        faceShareWithHex = .false.
        ! check settings are correct
    
        ! test if match with side 1
        if ((f1 == c2 .and. f2 == c1 .and. f3 == c4 .and. f4 == c3) .or. &
            (f1 == c1 .and. f2 == c4 .and. f3 == c3 .and. f4 == c2) .or. &
            (f1 == c4 .and. f2 == c3 .and. f3 == c2 .and. f4 == c1) .or. &
            (f1 == c3 .and. f2 == c2 .and. f3 == c1 .and. f4 == c4)) then
            faceShareWithHex = .true.
            return
        ! test if match with side 2
        elseif ((f1 == c6 .and. f2 == c5 .and. f3 == c1 .and. f4 == c2) .or. &
            (f1 == c5 .and. f2 == c1 .and. f3 == c2 .and. f4 == c6) .or. &
            (f1 == c1 .and. f2 == c2 .and. f3 == c6 .and. f4 == c5) .or. &
            (f1 == c2 .and. f2 == c6 .and. f3 == c5 .and. f4 == c1)) then
            faceShareWithHex = .true.
            return
        ! test if match with side 3
        elseif ((f1 == c7 .and. f2 == c6 .and. f3 == c2 .and. f4 == c3) .or. &
            (f1 == c3 .and. f2 == c7 .and. f3 == c6 .and. f4 == c2) .or. &
            (f1 == c2 .and. f2 == c3 .and. f3 == c7 .and. f4 == c6) .or. &
            (f1 == c6 .and. f2 == c2 .and. f3 == c3 .and. f4 == c7)) then
            faceShareWithHex = .true.
            return
        ! test if match with side 4
        elseif ((f1 == c8 .and. f2 == c7 .and. f3 == c3 .and. f4 == c4) .or. &
            (f1 == c7 .and. f2 == c3 .and. f3 == c4 .and. f4 == c8) .or. &
            (f1 == c3 .and. f2 == c4 .and. f3 == c8 .and. f4 == c7) .or. &
            (f1 == c4 .and. f2 == c8 .and. f3 == c7 .and. f4 == c3)) then
            faceShareWithHex = .true.
            return
        ! test if match with side 5
        elseif ((f1 == c5 .and. f2 == c8 .and. f3 == c4 .and. f4 == c1) .or. &
            (f1 == c8 .and. f2 == c4 .and. f3 == c1 .and. f4 == c5) .or. &
            (f1 == c4 .and. f2 == c1 .and. f3 == c5 .and. f4 == c8) .or. &
            (f1 == c1 .and. f2 == c5 .and. f3 == c8 .and. f4 == c4)) then
            faceShareWithHex = .true.
            return
        ! test if match with side 6
        elseif ((f1 == c5 .and. f2 == c6 .and. f3 == c7 .and. f4 == c8) .or. &
            (f1 == c6 .and. f2 == c7 .and. f3 == c8 .and. f4 == c5) .or. &
            (f1 == c7 .and. f2 == c8 .and. f3 == c5 .and. f4 == c6) .or. &
            (f1 == c8 .and. f2 == c5 .and. f3 == c6 .and. f4 == c7)) then
            faceShareWithHex = .true.
            return
        end if
    
    end function faceShareWithHex

    subroutine quadNormal(x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4,face_nrml_out,face_nrml_mag_out)
        use module_common_data, only : p2, half ! double precision
        implicit none
        real(p2), intent(in) :: x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4 ! vertices of quad face
        real(p2), dimension(3), intent(out) ::  face_nrml_out          ! unit normal vecotr
        real(p2),               intent(out) ::  face_nrml_mag_out      ! magnitude of normal vector (face area)

        ! Computed using the following link https://twitter.com/HiroNishikawa/status/1446403228258734081
        ! Twitter is obiously the best source of scholarly data
        ! Essentially the normal vector of a quad face is the sum of the two tri faces that combine to form it

        ! Some local vars
        ! real(p2), dimension(3) :: normal1, normal2      ! normal vectors of the two triangles
        ! real(p2)               :: area1, area2, mag     ! area of the two triangles
        ! real(p2)               :: nx, nx1
        
        ! call trinormal(x1,y1,z1, x2,y2,z2, x3,y3,z3, normal1, area1)
        ! call trinormal(x1,y1,z1, x3,y3,z3, x4,y4,z4, normal2, area2)
        ! nx1 = half*(y1*(z2-z3) + y2*(z3-z1) + y3*(z1-z2))
        ! face_nrml_out = normal1 + normal2
        ! mag = sqrt(face_nrml_out(1)**2 + face_nrml_out(2)**2 + face_nrml_out(3)**2)
        ! face_nrml_out = face_nrml_out/mag
        ! face_nrml_mag_out = area1 + area2
        face_nrml_out(1) = half * (y1*(z2-z4) + y2*(z3-z1) + y3*(z4-z2) + y4*(z1-z3))
        face_nrml_out(2) = half * (z1*(x2-x4) + z2*(x3-x1) + z3*(x4-x2) + z4*(x1-x3))
        face_nrml_out(3) = half * (x1*(y2-y4) + x2*(y3-y1) + x3*(y4-y2) + x4*(y1-y3))
        face_nrml_mag_out = sqrt(face_nrml_out(1)**2 + face_nrml_out(2)**2 + face_nrml_out(3)**2)
        face_nrml_out = face_nrml_out/face_nrml_mag_out
        ! write(*,*) face_nrml_out(1) - nx
        ! write(*,*) normal1(1) - nx1
        ! write(*,*)
    end subroutine quadNormal
    subroutine triNormal(x1,y1,z1, x2,y2,z2, x3,y3,z3,  normal, area)
        use module_common_data, only : p2, half ! double precision
        implicit none
        real(p2), intent(in) :: x1,y1,z1, x2,y2,z2, x3,y3,z3
        real(p2), dimension(3), intent(out) ::  normal          ! unit normal vecotr
        real(p2),               intent(out) ::  area      ! magnitude of normal vector (face area)
        ! Local Vars
        real(p2), dimension(3)              :: a, b 
        real(p2)                            :: mag
        real(p2) :: nx,ny,nz
        ! This should be cleaned up at some point, perhaps when I get to the chapter in my book on libraries...
        a(1) = x2-x1
        a(2) = y2-y1  
        a(3) = z2-z1
        b(1) = x3-x1
        b(2) = y3-y1  
        b(3) = z3-z1          
        ! ! Cross Product
        normal(1) = half*(a(2)*b(3) - a(3)*b(2))
        normal(2) = half*(a(3)*b(1) - a(1)*b(3))
        normal(3) = half*(a(1)*b(2) - a(2)*b(1))
        ! I guess Were using the above, it seems to have slightly better finite precision... :)

        ! normal(1) = half*(y1*(z2-z3)+y2*(z3-z1)+y3*(z1-z2))
        ! normal(2) = half*(z1*(x2-x3)+z2*(x3-x1)+z3*(x1-x2))
        ! normal(3) = half*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))
        ! write(*,*) normal(1) - nx
        ! write(*,*) normal(2) - ny
        ! write(*,*) normal(3) - nz
        area = sqrt(normal(1)**2 + normal(2)**2 + normal(3)**2)
        normal = normal / area
        ! area = half * mag        
    end subroutine triNormal
    subroutine get_tri_face_centroid(x1,y1,z1, x2,y2,z2, x3,y3,z3, triCentroid)
        use module_common_data, only : p2, one, three ! double precision
        implicit none
        real(p2), intent(in) :: x1,y1,z1, x2,y2,z2, x3,y3,z3  ! Three vertices describing a triangle
        real(p2), dimension(3), intent(out) ::  triCentroid   ! vector of cell centroid

        triCentroid(1) = (one/three)*(x1 + x2 + x3)
        triCentroid(2) = (one/three)*(y1 + y2 + y3)
        triCentroid(3) = (one/three)*(z1 + z2 + z3)
    end subroutine get_tri_face_centroid
    subroutine get_quad_face_centroid(x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4,  quadCentroid)
        use module_common_data, only : p2, one, four ! double precision
        implicit none
        real(p2), intent(in) :: x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4
        real(p2), dimension(3), intent(out) ::  quadCentroid
        
        quadCentroid(1) = (one/four)*(x1 + x2 + x3 + x4)
        quadCentroid(2) = (one/four)*(y1 + y2 + y3 + y4)
        quadCentroid(3) = (one/four)*(z1 + z2 + z3 + z4)
    end subroutine get_quad_face_centroid
    logical function checkInsideTet(x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4, xc,yc,zc)
        use module_common_data, only : p2           ! double precision
        implicit none
        real(p2), intent(in) :: x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4
        real(p2), intent(in) :: xc,yc,zc
        checkInsideTet = (samesides(x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4, xc,yc,zc) .and. &
                          samesides(x2,y2,z2, x3,y3,z3, x4,y4,z4, x1,y1,z1, xc,yc,zc) .and. &
                          samesides(x3,y3,z3, x4,y4,z4, x1,y1,z1, x2,y2,z2, xc,yc,zc) .and. &
                          samesides(x4,y4,z4, x1,y1,z1, x2,y2,z2, x3,y3,z3, xc,yc,zc))
    end function checkInsideTet
    logical function samesides(x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4, xc,yc,zc)
        use module_common_data, only : p2, one          ! double precision
        implicit none
        real(p2), intent(in) :: x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4
        real(p2), intent(in) :: xc,yc,zc

        ! Local Vars
        real(p2), dimension(3)              :: a, b, v4v1, P, normal 
        real(p2)                            :: dot4v, dotP
        ! This should be cleaned up at some point, perhaps when I get to the chapter in my book on libraries...
        a(1) = x2-x1
        a(2) = y2-y1  
        a(3) = z2-z1
        b(1) = x3-x1
        b(2) = y3-y1  
        b(3) = z3-z1          
        ! Cross Product
        normal(1) = (a(2)*b(3) - a(3)*b(2))
        normal(2) = (a(3)*b(1) - a(1)*b(3))
        normal(3) = (a(1)*b(2) - a(2)*b(1))
        v4v1 = (/x4-x1, y4-y1, z4-z1/)
        P = (/xc-x1, yc-y1, zc-z1/)
        dot4v = dot_product(normal,v4v1)
        dotP  = dot_product(normal,   P)
        samesides = ((dot4v < 0) .eqv. (dotP < 0)) ! check if both have the same sign, and don't have to worry about 
        ! overflow (x*y >= 0 could have issues for large x and y)
        ! note: treats 0 as positive
    end function samesides

    ! subroutine destroy_cells(ncells_to_be_destroyed,cell_to_be_destroyed)
    !     ! A quick subroutine to deallocate the components of a cell type that is no longer needed
    !     implicit none
    !     integer, intent(in) :: ncells_to_be_destroyed
    !     type(cc_data_type), intent(inout) :: cell_to_be_destroyed
    !     integer :: i

    !     do i = 1,ncells
    !         if (allocated(cell_to_be_destroyed(i)%commonv)) deallocate(cell_to_be_destroyed(i)%commonv)
    !         if (allocated(cell_to_be_destroyed(i)%vtx)) deallocate(cell_to_be_destroyed(i)%vtx)
    !         if (allocated(cell_to_be_destroyed(i)%nghbr)) deallocate(cell_to_be_destroyed(i)%nghbr)
    !         if (allocated(cell_to_be_destroyed(i)%ncommonv)) deallocate(cell_to_be_destroyed(i)%ncommonv)
    !         if (allocated(cell_to_be_destroyed(i)%gnghbr)) deallocate(cell_to_be_destroyed(i)%gnghbr)
            
            

    !     end do
    !     deallocate(cell_to_be_destroyed)
    ! end subroutine destroy_cells
end module module_ccfv_data_grid


