module module_ccfv_gradient
    
    use module_common_data, only : p2, one, zero
    implicit none
    
    ! Public data:
    public :: compute_lsq_coefficients
    public :: construct_vertex_stencil
    public :: compute_gradient
    public :: compute_temperature_gradient
    public :: my_alloc_int_ptr
    public :: perturbation_gradient
    public :: perturbation_gradient_boundary
    public :: compute_weighted_gradient

    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    ! Below are the definition of the data used to a LSQ gradient method.
    !------------------------------------------------------------------------------------

    !------------------------------------------
    !>> Cell-centered LSQ stencil data
    !------------------------------------------

    !Custom data type for lsq gradient stencil (nghbrs used in LSQ).
    type cc_lsq_data_type
        integer                           :: nnghbrs_lsq !number of lsq neighbors
        integer , dimension(:)  , pointer ::  nghbr_lsq  !list of lsq neighbors
        real(p2), dimension(:)  , pointer ::         cx  !LSQ coefficient for x-derivative
        real(p2), dimension(:)  , pointer ::         cy  !LSQ coefficient for y-derivative
        real(p2), dimension(:)  , pointer ::         cZ  !LSQ coefficient for z-derivative
    end type cc_lsq_data_type

    !Cell data array in the custom data type.
    type(cc_lsq_data_type), dimension(:), pointer :: cclsq  !cell-centered LSQ array
    type(cc_lsq_data_type), dimension(:), pointer :: cclsq_w  !cell-centered LSQ array

    !E.g., # of nghbrs at cell i = cclsq(i)%nnghbrs_lsq 
    !      a list of nghbrs at cell i = cclsq(i)%nghbr_lsq(1:cclsq(i)%nnghbrs_lsq)

    !These data will be allocated for a given grid size, and filled in the
    !following subroutine: compute_lsq_coefficients.

    !------------------------------------------------------------------------------------
    ! End of the data used to implement a LSQ gradient method.
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------

contains
    

    subroutine compute_gradient
        use module_common_data   , only : zero, p2
        use module_ccfv_data_grid, only : ncells
        use module_ccfv_data_soln, only : gradw, w, p_inf, gradq, q
        use module_input_parameter, only: gauge_pressure, low_mach_correction
        
        implicit none
        !real(p2), dimension(3,5,32)   :: local_gradw ! debugging
        real(p2) :: wi, wk, qi, qk
        integer  :: ivar, i, k, nghbr_cell
        
        if (low_mach_correction) then
            cell_loop_prim : do i = 1,ncells
                var_loop_prim : do ivar = 1,5
                    qi = q(ivar,i)
                    nghbr_loop_prim : do k = 1,cclsq(i)%nnghbrs_lsq
                        nghbr_cell = cclsq(i)%nghbr_lsq(k)
                        qk = q(ivar,nghbr_cell)
                        gradq(1,ivar,i) = gradq(1,ivar,i) + cclsq(i)%cx(k)*(qk-qi)
                        gradq(2,ivar,i) = gradq(2,ivar,i) + cclsq(i)%cy(k)*(qk-qi)
                        gradq(3,ivar,i) = gradq(3,ivar,i) + cclsq(i)%cz(k)*(qk-qi)
                    end do nghbr_loop_prim
                end do var_loop_prim
            end do cell_loop_prim
        else
            cell_loop : do i = 1,ncells
                var_loop : do ivar = 1,5
                    wi = w(ivar,i)
                    if (ivar == 5) wi = wi - p_inf + gauge_pressure
                    nghbr_loop : do k = 1,cclsq(i)%nnghbrs_lsq
                        nghbr_cell = cclsq(i)%nghbr_lsq(k)
                        wk = w(ivar,nghbr_cell)
                        if (ivar == 5) wk = wk - p_inf + gauge_pressure
                        gradw(1,ivar,i) = gradw(1,ivar,i) + cclsq(i)%cx(k)*(wk-wi)
                        gradw(2,ivar,i) = gradw(2,ivar,i) + cclsq(i)%cy(k)*(wk-wi)
                        gradw(3,ivar,i) = gradw(3,ivar,i) + cclsq(i)%cz(k)*(wk-wi)
                    end do nghbr_loop
                end do var_loop
            end do cell_loop
        end if

        ! compute the gradient for each of the primitive variables
    end subroutine compute_gradient
    
    subroutine compute_weighted_gradient
        use module_common_data   , only : zero, p2
        use module_ccfv_data_grid, only : ncells
        use module_ccfv_data_soln, only : gradw_w, w, p_inf, gradq_w, q
        use module_input_parameter, only: gauge_pressure, low_mach_correction
        
        implicit none
        !real(p2), dimension(3,5,32)   :: local_gradw ! debugging
        real(p2) :: wi, wk, qi, qk
        integer  :: ivar, i, k, nghbr_cell
        
        if (low_mach_correction) then
            cell_loop_prim : do i = 1,ncells
                var_loop_prim : do ivar = 1,5
                    qi = q(ivar,i)
                    nghbr_loop_prim : do k = 1,cclsq(i)%nnghbrs_lsq
                        nghbr_cell = cclsq(i)%nghbr_lsq(k)
                        qk = q(ivar,nghbr_cell)
                        gradq_w(1,ivar,i) = gradq_w(1,ivar,i) + cclsq(i)%cx(k)*(qk-qi)
                        gradq_w(2,ivar,i) = gradq_w(2,ivar,i) + cclsq(i)%cy(k)*(qk-qi)
                        gradq_w(3,ivar,i) = gradq_w(3,ivar,i) + cclsq(i)%cz(k)*(qk-qi)
                    end do nghbr_loop_prim
                end do var_loop_prim
            end do cell_loop_prim
        else
            cell_loop : do i = 1,ncells
                var_loop : do ivar = 1,5
                    wi = w(ivar,i)
                    if (ivar == 5) wi = wi - p_inf + gauge_pressure
                    nghbr_loop : do k = 1,cclsq(i)%nnghbrs_lsq
                        nghbr_cell = cclsq(i)%nghbr_lsq(k)
                        wk = w(ivar,nghbr_cell)
                        if (ivar == 5) wk = wk - p_inf + gauge_pressure
                        gradw_w(1,ivar,i) = gradw_w(1,ivar,i) + cclsq(i)%cx(k)*(wk-wi)
                        gradw_w(2,ivar,i) = gradw_w(2,ivar,i) + cclsq(i)%cy(k)*(wk-wi)
                        gradw_w(3,ivar,i) = gradw_w(3,ivar,i) + cclsq(i)%cz(k)*(wk-wi)
                    end do nghbr_loop
                end do var_loop
            end do cell_loop
        end if

        ! compute the gradient for each of the primitive variables
    end subroutine compute_weighted_gradient

    subroutine perturbation_gradient(gradw_orig1,gradw_orig2,ivar,ci,ck,u_perturb,grad_perturb1,grad_perturb2)
        use module_common_data , only : zero, p2
        use module_ccfv_data_soln, only : w, u2w

        implicit none
        real(p2), dimension(3,5), intent(in)    :: gradw_orig1,gradw_orig2          ! Original input gradient
        integer                 , intent(in)    :: ivar                             ! Variable that is being perturbed
        integer,                  intent(in)    :: ci,ck                            ! Cell that is being perturbed
        real(p2), dimension(5),   intent(in)    :: u_perturb
        real(p2), dimension(3,5), intent(out)   :: grad_perturb1,grad_perturb2      ! Resultant perturbed gradient

        ! Local vars
        real(p2), dimension(5) :: w_perturb
        real(p2) :: wi, wk
        integer  :: i,k,nghbr_cell

        grad_perturb1 = gradw_orig1
        grad_perturb1(:,ivar) = zero
        grad_perturb2 = gradw_orig2
        grad_perturb2(:,ivar) = zero

        w_perturb = u2w(u_perturb)
        ! Update gradient at cell i
        wi = w_perturb(ivar)
        nghbr_loop : do k = 1,cclsq(ci)%nnghbrs_lsq
            nghbr_cell = cclsq(ci)%nghbr_lsq(k)
            wk = w(ivar,nghbr_cell)
            grad_perturb1(1,ivar) = grad_perturb1(1,ivar) + cclsq(ci)%cx(k)*(wk-wi)
            grad_perturb1(2,ivar) = grad_perturb1(2,ivar) + cclsq(ci)%cy(k)*(wk-wi)
            grad_perturb1(3,ivar) = grad_perturb1(3,ivar) + cclsq(ci)%cz(k)*(wk-wi)
            !local_gradw = gradw
        end do nghbr_loop

        ! update gradient at cell k
        wi = w(ivar,ck)
        nghbr_loopk : do k = 1,cclsq(ck)%nnghbrs_lsq
            nghbr_cell = cclsq(ck)%nghbr_lsq(k)
            if (nghbr_cell == ci) then
                wk = w_perturb(ivar)
            else
                wk = w(ivar,nghbr_cell)
            end if
            grad_perturb2(1,ivar) = grad_perturb2(1,ivar) + cclsq(ck)%cx(k)*(wk-wi)
            grad_perturb2(2,ivar) = grad_perturb2(2,ivar) + cclsq(ck)%cy(k)*(wk-wi)
            grad_perturb2(3,ivar) = grad_perturb2(3,ivar) + cclsq(ck)%cz(k)*(wk-wi)
            !local_gradw = gradw
        end do nghbr_loopk

    end subroutine perturbation_gradient

    subroutine perturbation_gradient_boundary(gradw_orig1,ivar,ci,u_perturb,grad_perturb1)
        use module_common_data , only : zero, p2
        use module_ccfv_data_soln, only : w, u2w

        implicit none
        real(p2), dimension(3,5), intent(in)    :: gradw_orig1          ! Original input gradient
        integer                 , intent(in)    :: ivar                             ! Variable that is being perturbed
        integer,                  intent(in)    :: ci                            ! Cell that is being perturbed
        real(p2), dimension(5),   intent(in)    :: u_perturb
        real(p2), dimension(3,5), intent(out)   :: grad_perturb1      ! Resultant perturbed gradient

        ! Local vars
        real(p2), dimension(5) :: w_perturb
        real(p2) :: wi, wk
        integer  :: i,k,nghbr_cell

        grad_perturb1 = gradw_orig1
        grad_perturb1(:,ivar) = zero


        w_perturb = u2w(u_perturb)
        ! Update gradient at cell i
        wi = w_perturb(ivar)
        nghbr_loop : do k = 1,cclsq(ci)%nnghbrs_lsq
            nghbr_cell = cclsq(ci)%nghbr_lsq(k)
            wk = w(ivar,nghbr_cell)
            grad_perturb1(1,ivar) = grad_perturb1(1,ivar) + cclsq(ci)%cx(k)*(wk-wi)
            grad_perturb1(2,ivar) = grad_perturb1(2,ivar) + cclsq(ci)%cy(k)*(wk-wi)
            grad_perturb1(3,ivar) = grad_perturb1(3,ivar) + cclsq(ci)%cz(k)*(wk-wi)

        end do nghbr_loop

        
    end subroutine perturbation_gradient_boundary

    subroutine perturb_temperature_gradient(gradT_orig1,gradT_orig2,ci,ck,T_perturb,grad_perturb1,grad_perturb2)
        use module_common_data, only : p2
        use module_ccfv_data_grid, only : ncells
        use module_ccfv_data_soln, only : w, Temp, gamma
        implicit none
        real(p2), dimension(3), intent(in)      :: gradT_orig1,gradT_orig2
        integer,                intent(in)      :: ci, ck
        real(p2),               intent(in)      :: T_perturb
        real(p2), dimension(3), INTENT(OUT)     :: grad_perturb1,grad_perturb2

        integer :: i, nghbr_cell, k
        real(p2) :: Ti, Tk
    
        ! Now we use the least squares stuff to calculate the gradient
        grad_perturb1 = zero
        
        nghbr_loop : do k = 1,cclsq(ci)%nnghbrs_lsq
            nghbr_cell = cclsq(ci)%nghbr_lsq(k)
            grad_perturb1(1) = grad_perturb1(1) + cclsq(i)%cx(k)*(Temp(k) - T_perturb)
            grad_perturb1(2) = grad_perturb1(2) + cclsq(i)%cy(k)*(Temp(k) - T_perturb)
            grad_perturb1(3) = grad_perturb1(3) + cclsq(i)%cz(k)*(Temp(k) - T_perturb)
        end do nghbr_loop
        
        ! update gradient at cell k
        Ti = Temp(ck)
        nghbr_loopk : do k = 1,cclsq(ck)%nnghbrs_lsq
            nghbr_cell = cclsq(ck)%nghbr_lsq(k)
            if (nghbr_cell == ci) then
                Tk = T_perturb
            else
                Tk = Temp(nghbr_cell)
            end if
            grad_perturb2(1) = grad_perturb2(1) + cclsq(ck)%cx(k)*(Tk-Ti)
            grad_perturb2(2) = grad_perturb2(2) + cclsq(ck)%cy(k)*(Tk-Ti)
            grad_perturb2(3) = grad_perturb2(3) + cclsq(ck)%cz(k)*(Tk-Ti)
            !local_gradw = gradw
        end do nghbr_loopk
    end subroutine perturb_temperature_gradient

    subroutine compute_temperature_gradient
        use module_common_data, only : p2
        use module_ccfv_data_grid, only : ncells
        use module_ccfv_data_soln, only : w, Temp, gamma, gradT
        
        implicit none
        integer :: i, nghbr_cell, k
    
        ! First we need to compute T
        do i = 1,ncells
            Temp(i) = w(5,i) * gamma / w(1,i)
        end do

        ! Now we use the least squares stuff to calculate the gradient
        gradT = zero
        cell_loop : do i = 1,ncells
            nghbr_loop : do k = 1,cclsq(i)%nnghbrs_lsq
                nghbr_cell = cclsq(i)%nghbr_lsq(k)
                gradT(1,i) = gradT(1,i) + cclsq(i)%cx(k)*(Temp(k) - Temp(i))
                gradT(2,i) = gradT(2,i) + cclsq(i)%cy(k)*(Temp(k) - Temp(i))
                gradT(3,i) = gradT(3,i) + cclsq(i)%cz(k)*(Temp(k) - Temp(i))
            end do nghbr_loop
        end do cell_loop
    end subroutine compute_temperature_gradient

    subroutine construct_vertex_stencil
        use module_common_data   , only : nnodes
        use module_ccfv_data_grid, only : cell, ncells
        use module_input_parameter,only : navier_stokes

        implicit none

        !To store the list of cells around each node.
        type node_type
            integer                        :: nc
            integer, dimension(:), pointer :: c
        end type node_type

        type(node_type), dimension(:), pointer :: node

        integer :: i, k, vk, kc, ii
        integer :: candidate_cell, n
        logical :: already_added

        integer :: ave_nghbr, min_nghbr, max_nghbr !<- For statistics
        integer :: max_nghbrs, imin, imax, icount  !<- For statistics
        integer, dimension(100) :: nghbrs !<- For statistics. Let's assume max # of nghbrs is 100. 

        write(*,*)
        write(*,*) " --------------------------------------------------"
        write(*,*) " Constructing vertex neighbors... "
        write(*,*)

        ! Initialization for statistical quantities.

        max_nghbrs = 0
        nghbrs = 0

        ave_nghbr = 0
        min_nghbr = 10000
        max_nghbr =-10000
             imin = 1
             imax = 1

        ! First, create node-to-cell lists for convenience.
        ! [Same as done in edu2d_module_ccfv_data_grid.f90.]

        allocate( node(nnodes) )

        do i = 1,nnodes
            node(i)%nc = 0
        end do

        do i = 1, ncells
            do k = 1, cell(i)%nvtx
                         vk = cell(i)%vtx(k)
                node(vk)%nc = node(vk)%nc + 1
            end do
        end do

        do i = 1,nnodes
            allocate( node(i)%c( node(i)%nc ) )
        end do
        do i = 1,nnodes
            node(i)%nc = 0
        end do
        do i = 1, ncells
            do k = 1, cell(i)%nvtx
                                     vk = cell(i)%vtx(k)
                node(vk)%nc             = node(vk)%nc + 1
                node(vk)%c(node(vk)%nc) = i
            end do
        end do

        !Allocate and initialize the LSQ (custom data) array.

        allocate(cclsq(ncells))
        if (navier_stokes) allocate(cclsq_w(ncells))

        ! Initialize it
        do i = 1,ncells
            cclsq(i)%nnghbrs_lsq = 0
        end do

        !------------------------------------------------------------------------
        ! Find vertex neighbors at each cell.

        cell_loop : do i = 1, ncells
            cell_vertex_loop : do k = 1, cell(i)%nvtx
                vk = cell(i)%vtx(k)
                node2cell_loop : do kc = 1, node(vk)%nc
                    candidate_cell = node(vk)%c(kc)
                    if (candidate_cell == i) cycle node2cell_loop !Skip the cell i.
                    ! Check if cell is already added
                    already_added = .false.
                    if (cclsq(i)%nnghbrs_lsq > 0) then
                        do ii = 1, cclsq(i)%nnghbrs_lsq
                            if ( candidate_cell == cclsq(i)%nghbr_lsq(ii) ) then
                                already_added = .true.
                                exit !OK, it is already added. Go to the next candidate.
                            endif
                        end do
                    endif

                    ! If candidate_cell is not yet added add it
                    if (.not.already_added) then
                                            n = cclsq(i)%nnghbrs_lsq + 1
                         cclsq(i)%nnghbrs_lsq = n                    !Increase the size by 1.
                        if (n==1) then
                            if (associated(cclsq(i)%nghbr_lsq)) nullify(cclsq(i)%nghbr_lsq)
                            allocate(cclsq(i)%nghbr_lsq(1))
                        else
                            call my_alloc_int_ptr(cclsq(i)%nghbr_lsq, n) !Expand the array by 1.
                        end if
                        cclsq(i)%nghbr_lsq(n) = candidate_cell       !Store the candidate.
                    end if
                end do node2cell_loop
            end do cell_vertex_loop
            ! Allocate the LSQ coefficient arrays for the cell(i)
            allocate( cclsq(i)%cx(cclsq(i)%nnghbrs_lsq) )
            allocate( cclsq(i)%cy(cclsq(i)%nnghbrs_lsq) )
            allocate( cclsq(i)%cz(cclsq(i)%nnghbrs_lsq) )

            ! record some statistics
            ave_nghbr = ave_nghbr + cclsq(i)%nnghbrs_lsq
            if (cclsq(i)%nnghbrs_lsq < min_nghbr) imin = i
            if (cclsq(i)%nnghbrs_lsq > max_nghbr) imax = i
            min_nghbr = min(min_nghbr, cclsq(i)%nnghbrs_lsq)
            max_nghbr = max(max_nghbr, cclsq(i)%nnghbrs_lsq)
        
            icount         = cclsq(i)%nnghbrs_lsq
            max_nghbrs     = max(max_nghbrs, icount)
            nghbrs(icount) = nghbrs(icount) + 1
        end do cell_loop

        ! Allocate weighted lsq structure if navier_stokes
        if (navier_stokes) cclsq_w = cclsq

        ! Deallocate 'node' as we don't need it any more.

        deallocate( node )
        
        write(*,*)
        write(*,*) " Constructed vertex neighbors "
        write(*,*)
        write(*,*) "      ave_nghbr = ", ave_nghbr/ncells
        write(*,*) "      min_nghbr = ", min_nghbr, " elm = ", imin
        write(*,*) "      max_nghbr = ", max_nghbr, " elm = ", imax
        write(*,*)
        do k = 1, max_nghbrs
            write(*,*) " # of neighbors = ", k,  ": # of such cells = ", nghbrs(k)
        end do
        write(*,*)
        write(*,*) " End of Constructing vertex neighbors... "
        write(*,*) " --------------------------------------------------"
        write(*,*)
    end subroutine construct_vertex_stencil

    subroutine compute_lsq_coefficients
        use module_common_data   , only : p2, zero, one, two
        use module_ccfv_data_grid, only : cell, ncells
        use module_input_parameter,only : navier_stokes

        implicit none

        integer                           :: i, k, nghbr_cell
        integer                           :: m, n             !Size of LSQ matrix: A(m,n).
        real(p2), pointer, dimension(:,:) :: a                !LSQ matrix: A(m,n).
        real(p2), pointer, dimension(:,:) :: rinvqt           !Pseudo inverse R^{-1}*Q^T
        real(p2)                          :: dx, dy, dz, weight_k!,lsq_weight_invdis_power
        real(p2)                          :: maxdx, maxdy, maxdz
        integer                           :: ix, iy, iz, lsq_weight_invdis_power
        real(p2), dimension(3)            :: maxDeltasNZ

        real(p2)                          :: xk, yk, zk, xi, yi, zi, wx, wy, wz
        logical                           :: verification_error
    
         !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        ! Allocate the custom-data LSQ array.

        call construct_vertex_stencil

        write(*,*)
        write(*,*) "--------------------------------------------------"
        write(*,*) " Computing LSQ coefficients... "
        write(*,*)
        
        ! Just for convenience.
    
        ix = 1
        iy = 2
        iz = 3
        maxdx = zero
        maxdy = zero
        maxdz = zero

        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        ! The power to the inverse distance weight. The value 0.0 is used to avoid
        ! instability known for Euler solvers. So, this is the unweighted LSQ gradient.
        ! More accurate gradients are obtained with 1.0, and such can be used for the
        ! viscous terms and source terms in turbulence models.
        lsq_weight_invdis_power = 0
         
        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        ! Compute the LSQ coefficients (cx,cy,cz) in all cells.
        cell_loop : do i = 1,ncells
            ! Define the LSQ problem size
            m = cclsq(i)%nnghbrs_lsq  ! # of nghbrs
            n = 3                     ! # of derivatives (unknowns = 3, (ux,uy,uz) in 3D)
            ! Allocate LSQ matrix and the pseudo inverse, R^{-1}*Q^T
            allocate(a(m,n)) ! note: it may produce some additional speed to switch the rows and columns here
            ! however a is a very small matrix and this subroutine is called once so we'll leave that for another day
            allocate(rinvqt(n,m))
            ! Initialize a
            a = zero
            !-------------------------------------------------------
            ! Build the weighted-LSQ matrix A(m,n).
            !
            !     weight_1 * [ (x1-xi)*wxi + (y1-yi)*wyi + (z1-zi)*wzi ] = weight_1 * [ w1 - wi ]
            !     weight_2 * [ (x2-xi)*wxi + (y2-yi)*wyi + (z2-zi)*wzi ] = weight_2 * [ w2 - wi ]
            !                 .
            !                 .
            !     weight_m * [ (xm-xi)*wxi + (ym-yi)*wyi + (zm-zi)*wzi ] = weight_2 * [ wm - wi ]
            nghbr_loop : do k = 1, m
                nghbr_cell = cclsq(i)%nghbr_lsq(k) !Neighbor cell number
                dx = cell(nghbr_cell)%xc - cell(i)%xc
                dy = cell(nghbr_cell)%yc - cell(i)%yc
                dz = cell(nghbr_cell)%zc - cell(i)%zc
                weight_k = one / sqrt( dx**2 + dy**2 + dz**2 )**lsq_weight_invdis_power
                a(k,1) = weight_k*dx
                a(k,2) = weight_k*dy
                a(k,3) = weight_k*dz
                maxdx  = max(abs(dx),maxdx)
                maxdy  = max(abs(dy),maxdy)
                maxdz  = max(abs(dz),maxdz)
            end do nghbr_loop
            !-------------------------------------------------------
            ! Perform QR factorization and compute R^{-1}*Q^T from A(m,n).
            call qr_factorization(a,rinvqt,m,n)

            !-------------------------------------------------------
            ! Compute and store the LSQ coefficients: R^{-1}*Q^T*w.
            !
            ! (wx,wy,wz) = R^{-1}*Q^T*RHS
            !            = sum_k (cx,cy,cz)*(wk-wi).

            nghbr_loop2 : do k = 1, m
                nghbr_cell = cclsq(i)%nghbr_lsq(k)
                dx = cell(nghbr_cell)%xc - cell(i)%xc
                dy = cell(nghbr_cell)%yc - cell(i)%yc
                dz = cell(nghbr_cell)%zc - cell(i)%zc
                weight_k = one / sqrt( dx**2 + dy**2 + dz**2 )**lsq_weight_invdis_power
                cclsq(i)%cx(k)  = rinvqt(ix,k) * weight_k
                cclsq(i)%cy(k)  = rinvqt(iy,k) * weight_k
                cclsq(i)%cz(k)  = rinvqt(iz,k) * weight_k
            end do nghbr_loop2
            !-------------------------------------------------------
            ! Deallocate a and rinvqt, whose size may change in the next cell. 
            deallocate(a, rinvqt)
        end do cell_loop

        lsq_weight_invdis_power = 1

        if (navier_stokes) then
            !--------------------------------------------------------------------------------
            !--------------------------------------------------------------------------------
            ! Compute the LSQ coefficients (cx,cy,cz) in all cells.
            cell_loop_w : do i = 1,ncells
                ! Define the LSQ problem size
                m = cclsq_w(i)%nnghbrs_lsq  ! # of nghbrs
                n = 3                     ! # of derivatives (unknowns = 3, (ux,uy,uz) in 3D)
                ! Allocate LSQ matrix and the pseudo inverse, R^{-1}*Q^T
                allocate(a(m,n)) ! note: it may produce some additional speed to switch the rows and columns here
                ! however a is a very small matrix and this subroutine is called once so we'll leave that for another day
                allocate(rinvqt(n,m))
                ! Initialize a
                a = zero
                !-------------------------------------------------------
                ! Build the weighted-LSQ matrix A(m,n).
                !
                !     weight_1 * [ (x1-xi)*wxi + (y1-yi)*wyi + (z1-zi)*wzi ] = weight_1 * [ w1 - wi ]
                !     weight_2 * [ (x2-xi)*wxi + (y2-yi)*wyi + (z2-zi)*wzi ] = weight_2 * [ w2 - wi ]
                !                 .
                !                 .
                !     weight_m * [ (xm-xi)*wxi + (ym-yi)*wyi + (zm-zi)*wzi ] = weight_2 * [ wm - wi ]
                nghbr_loop_w : do k = 1, m
                    nghbr_cell = cclsq_w(i)%nghbr_lsq(k) !Neighbor cell number
                    dx = cell(nghbr_cell)%xc - cell(i)%xc
                    dy = cell(nghbr_cell)%yc - cell(i)%yc
                    dz = cell(nghbr_cell)%zc - cell(i)%zc
                    weight_k = one / sqrt( dx**2 + dy**2 + dz**2 )**lsq_weight_invdis_power
                    a(k,1) = weight_k*dx
                    a(k,2) = weight_k*dy
                    a(k,3) = weight_k*dz
                    maxdx  = max(abs(dx),maxdx)
                    maxdy  = max(abs(dy),maxdy)
                    maxdz  = max(abs(dz),maxdz)
                end do nghbr_loop_w
                !-------------------------------------------------------
                ! Perform QR factorization and compute R^{-1}*Q^T from A(m,n).
                call qr_factorization(a,rinvqt,m,n)

                !-------------------------------------------------------
                ! Compute and store the LSQ coefficients: R^{-1}*Q^T*w.
                !
                ! (wx,wy,wz) = R^{-1}*Q^T*RHS
                !            = sum_k (cx,cy,cz)*(wk-wi).

                nghbr_loop2_w : do k = 1, m
                    nghbr_cell = cclsq_w(i)%nghbr_lsq(k)
                    dx = cell(nghbr_cell)%xc - cell(i)%xc
                    dy = cell(nghbr_cell)%yc - cell(i)%yc
                    dz = cell(nghbr_cell)%zc - cell(i)%zc
                    weight_k = one / sqrt( dx**2 + dy**2 + dz**2 )**lsq_weight_invdis_power
                    cclsq_w(i)%cx(k)  = rinvqt(ix,k) * weight_k
                    cclsq_w(i)%cy(k)  = rinvqt(iy,k) * weight_k
                    cclsq_w(i)%cz(k)  = rinvqt(iz,k) * weight_k
                end do nghbr_loop2_w
                !-------------------------------------------------------
                ! Deallocate a and rinvqt, whose size may change in the next cell. 
                deallocate(a, rinvqt)
            end do cell_loop_w
        endif




        ! Verification
        ! Compute the gradient of w = 2*x+y+4*z to se if we get wx = 2, wy = 1, and wz = 4 correctly
        verification_error = .false.

        do i = 1,ncells
            ! initialize wx, wy, and wz
            wx = zero
            wy = zero
            wz = zero
            ! (xi,yi,zi) to be used to compute the function 2*x+y+4z at i
            xi = cell(i)%xc
            yi = cell(i)%yc
            zi = cell(i)%zc

            ! look over the vertex neighboes
            do k = 1,cclsq(i)%nnghbrs_lsq
                nghbr_cell = cclsq(i)%nghbr_lsq(k)
                xk = cell(nghbr_cell)%xc
                yk = cell(nghbr_cell)%yc
                zk = cell(nghbr_cell)%zc
                ! This is how we use the LSQ coefficients: accumulate cx*(wk-wi)
                ! and cy*(wk-wi) and cz*(wk-wi)
                wx = wx + cclsq(i)%cx(k)*( (2.0*xk+yk+4.0*zk)-(2.0*xi+yi+4.0*zi) )
                wy = wy + cclsq(i)%cy(k)*( (2.0*xk+yk+4.0*zk)-(2.0*xi+yi+4.0*zi) )
                wz = wz + cclsq(i)%cz(k)*( (2.0*xk+yk+4.0*zk)-(2.0*xi+yi+4.0*zi) )
            end do
            maxDeltasNZ = zero
            if (maxdx > 0.001_p2) maxDeltasNZ(1) = one
            if (maxdy > 0.001_p2) maxDeltasNZ(2) = one
            if (maxdz > 0.001_p2) maxDeltasNZ(3) = one
            if ( maxDeltasNZ(1)*abs(wx-two) > 1.0e-06_p2 .or. &
                 maxDeltasNZ(2)*abs(wy-one) > 1.0e-06_p2 .or. &
                 maxDeltasNZ(3)*abs(wz-4.0_p2) > 1.0e-06_p2) then
                    write(*,*) " wx = ", wx, " exact ux = 2.0"!,maxDeltasNZ(1)*abs(wx-two)
                    write(*,*) " wy = ", wy, " exact uy = 1.0"!,maxDeltasNZ(2)*abs(wy-one)
                    write(*,*) " wz = ", wz, " exact uz = 4.0"!, maxDeltasNZ(3)*abs(wz-4.0_p2),maxDeltasNZ(3)
                    verification_error = .true.
            end if

            if (navier_stokes) then
                wx = zero
                wy = zero
                wz = zero
                ! look over the vertex neighboes
                do k = 1,cclsq_w(i)%nnghbrs_lsq
                    nghbr_cell = cclsq_w(i)%nghbr_lsq(k)
                    xk = cell(nghbr_cell)%xc
                    yk = cell(nghbr_cell)%yc
                    zk = cell(nghbr_cell)%zc
                    ! This is how we use the LSQ coefficients: accumulate cx*(wk-wi)
                    ! and cy*(wk-wi) and cz*(wk-wi)
                    wx = wx + cclsq_w(i)%cx(k)*( (2.0*xk+yk+4.0*zk)-(2.0*xi+yi+4.0*zi) )
                    wy = wy + cclsq_w(i)%cy(k)*( (2.0*xk+yk+4.0*zk)-(2.0*xi+yi+4.0*zi) )
                    wz = wz + cclsq_w(i)%cz(k)*( (2.0*xk+yk+4.0*zk)-(2.0*xi+yi+4.0*zi) )
                end do
                maxDeltasNZ = zero
                if (maxdx > 0.001_p2) maxDeltasNZ(1) = one
                if (maxdy > 0.001_p2) maxDeltasNZ(2) = one
                if (maxdz > 0.001_p2) maxDeltasNZ(3) = one
                if ( maxDeltasNZ(1)*abs(wx-two) > 1.0e-06_p2 .or. &
                    maxDeltasNZ(2)*abs(wy-one) > 1.0e-06_p2 .or. &
                    maxDeltasNZ(3)*abs(wz-4.0_p2) > 1.0e-06_p2) then
                        write(*,*) " weighted wx = ", wx, " exact ux = 2.0"!,maxDeltasNZ(1)*abs(wx-two)
                        write(*,*) " weighted wy = ", wy, " exact uy = 1.0"!,maxDeltasNZ(2)*abs(wy-one)
                        write(*,*) " weighted wz = ", wz, " exact uz = 4.0"!, maxDeltasNZ(3)*abs(wz-4.0_p2),maxDeltasNZ(3)
                        verification_error = .true.
                end if
            endif
        end do
        if (verification_error) then

            write(*,*) " LSQ coefficients are not correct. See above. Stop."
            stop
         
        else
         
            write(*,*) " Verified: LSQ coefficients are exact for a linear function."
         
        endif
         
        write(*,*)
        write(*,*) " End of Computing LSQ coefficients... "
        write(*,*) "--------------------------------------------------"
        write(*,*)
        
        ! End of Compute the LSQ coefficients in all cells.
        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
         
    end subroutine compute_lsq_coefficients
    !****************************************************************************
    ! ------------------ QR Factorization ---------------------
    !
    !  This subroutine solves the LSQ problem: A*x=b, A=mxn matrix.
    !
    !  IN :       a = (m x n) LSQ matrix.  (m >= n)
    !
    ! OUT :  rinvqt = R^{-1}*Q^t, which gives the solution as x = R^{-1}*Q^t*b.
    !
    !*****************************************************************************
    subroutine qr_factorization(a,rinvqt,m,n)

        implicit none
    
        integer , parameter ::   p2 = selected_real_kind(P=15) !Double precision
        real(p2), parameter :: zero = 0.0_p2
        real(p2), parameter ::  one = 1.0_p2
        real(p2), parameter ::  two = 2.0_p2
    
        !Input
        integer ,                 intent( in) :: m, n
        real(p2), dimension(m,n), intent( in) :: a
    
        !Output
        real(p2), dimension(n,m), intent(out) :: rinvqt
            
        ! Local variables
        !
        ! Note: Think if you can reduce the number of
        !       variables below to save memory.
        
        integer                  :: i, j, k, ii
    
        real(p2), dimension(:,:), pointer :: r
        real(p2)                          :: abs_rk, sign_rkk, wTw
        real(p2), dimension(m)            :: w, rk
        real(p2), dimension(m,m)          :: qt, wwT
    
        real(p2), dimension(:,:), pointer ::  r_nxn
        real(p2), dimension(:),   pointer   ::   y, b
        real(p2)                 ::    rhs
        integer,  dimension(n)   :: zeroCheck
        integer                  :: nonZeroColumns
        real(p2), dimension(n,m) :: rinvqt_intermediate
        
        do i = 1,n
            zeroCheck(i) = 0
        end do
        nonZeroColumns = 0
        do i = 1,n
            if (.not.(vector_norm(a(:,i),m)) == 0) then
                zeroCheck(i) = 1
                nonZeroColumns = nonZeroColumns + 1
            end if
        end do
        if (nonZeroColumns == 0) then
            write(*,*) " error dx = dy = dz = 0" ! This really shouldn't happen...
        end if
        allocate(r(m,nonZeroColumns))
        ii = 1
        do i = 1,n
            if (zeroCheck(i) == 1) then
                r(:,ii) = a(:,i)
                ii = ii + 1
            end if
        end do

        if (m < n) then
        write(*,*) " Underdetermined system detected... m < n: "
        write(*,*) "   m =  ", m
        write(*,*) "   n =  ", n
        write(*,*) " qr_factorization() not designed to solve such a problem... Stop. "
        stop
        endif
    
        !-------------------------------------------------------
        ! Initialization: R = A
    
        ! r = a ! not anymore
    
        !-------------------------------------------------------
        ! Initialization: Qt = I
        
        qt = zero
    
        do i = 1, m
            qt(i,i) = one
        end do
    
        !-------------------------------------------------------
        ! Apply reflection to each column of R, and generate
        ! the final upper triangular matrix R and the transpose
        ! Qt of the orthonormal matrix Q.
    
        column_loop : do k = 1, nonZeroColumns
    
            !Our target are the elements below the (k,k) element
            !in k-th column, i.e., r(k:m).
            !So, rk gets shorter as we move on (as k increases).
    
            rk      = zero
            rk(k:m) = r(k:m,k)
    
            !Reflector Hk will zero out all the elements below r(k).
        
            !Compute the length of rk and the sign of the kth element.
    
              abs_rk = sqrt( dot_product(rk,rk) )
            sign_rkk = sign( one, rk(k) )
    
            !Define the reflecting vector w:   w = |rk|*(1,0,0,...,0)-rk
            !                               or w =-|rk|*(1,0,0,...,0)-rk
            !We switch the reflection (there are two possible ones)
            !to avoid w = 0 that can happen if rk=(1,0,0,...,0).
    
            w      = zero
            w(k)   = -sign_rkk*abs_rk
            w(k:m) = w(k:m) - rk(k:m)
    
            !Compute the length^2 of w: wt*w = [x,x,...,x]|x| = dot product = scalar.
            !                                             |x|
            !                                             |.|
            !                                             |.|
            !                                             |x|
        
            wTw = dot_product(w,w)
        
            !Compute the dyad of w: w*wt = |x|[x,x,...,x] = mxm matrix.
            !                              |x|
            !                              |.|
            !                              |.|
            !                              |x|
        
            do i = 1, m
                do j = 1, m
                    wwT(i,j) = w(i)*w(j)
                end do
            end do
        
            !We now apply the reflector matrix Hk = I-2*wwt/wTw,
            !and update R and Qt.
        
            !Update  R:  R = Hk*R  = (I-2*wwt/wTw)*R  = R-2*(wwt*R)/wTw
        
            r  =  r - two*matmul(wwT,r)/wTw
        
            !Update Qt: Qt = Hk*Qt = (I-2*wwt/wTw)*Qt = Qt-2*(wwt*Qt)/wTw
        
            qt = qt - two*matmul(wwT,qt)/wTw
        
        end do column_loop
    
        !-------------------------------------------------------
        ! Compute rinvqt(1:n,1:m) = R_{nxn}^{-1} * Q_{nxm}^t by
        ! solving R_{nxn} * rinvqt(1:n,k) = Q_{nxm}^t(1:n,k)
        ! for k=1,n. We can solve it easily by back substitution
        ! since R_{nxn} is upper triangular.
        
        allocate(r_nxn(nonZeroColumns,nonZeroColumns))
        r_nxn =  r(1:nonZeroColumns,1:nonZeroColumns)
    
        allocate(y(nonZeroColumns))
        allocate(b(nonZeroColumns))
        rinvqt = zero ! initialize as zero since some values won't get assigned later
        rinvqt_intermediate = zero
        do k = 1, m
    
            !Solve r*y = b, where y is the k-th column of rinvqt.
    
            b = qt(1:nonZeroColumns,k)
    
            !Solve the lower right equation.
    
            rhs = b(nonZeroColumns)
            y(nonZeroColumns) = rhs/r_nxn(nonZeroColumns,nonZeroColumns)
    
            !Go up and solve.
            do i = nonZeroColumns-1, 1, -1
    
                !Take all known parts (j=i+1,n) to the rhs.
        
                !RHS is known, of course.
                rhs = b(i)
                !Below are all known since the solutions y(j=i+1,n) has already been solved.
                do j = i+1, nonZeroColumns
                    rhs = rhs - r_nxn(i,j)*y(j)
                end do
        
                !Divide the rhs by the coefficient of the (i,i) part.
                y(i) = rhs/r_nxn(i,i)
        
            end do
    
            !The soluton x is the k-th column of rinvqt.
            rinvqt_intermediate(:nonZeroColumns,k) = y(:)    
        end do
        nonZeroColumns = 0
        do i = 1,n
            if (zeroCheck(i) == 1) then
                nonZeroColumns = nonZeroColumns + 1
                rinvqt(i,:) = rinvqt_intermediate(nonZeroColumns,:)
            end if
        end do   
        deallocate(r,r_nxn,y,b) 
    end subroutine qr_factorization
    ! Function for computing the L2 vector norm from http://me.rice.edu/~akin/OOP_Copyrighted/4_Features_of_Lang/vector_norm.f90
    ! modified from a standalone program
    real function vector_norm(x, n)
        !   A simple vector norm program, Fortran 90 version
        use module_common_data   , only : p2, zero, one, two  
        implicit none
            !integer, parameter :: dp = selected_real_kind(14) ! 14 digits
            real(p2), dimension(:), intent(in)  :: x    
            integer, intent(in)                 :: n
            !integer, parameter :: MAX = 1000;           !  Maximum vector size
            !integer n;                                  !  Actual vector size
            !integer i;                                  !  Loop controls
            !real(p2) :: x(MAX)                          !  The vector
            !real(p2) :: L1_norm, L2_norm;               !  Work results
            
            !print *,"Enter the vector length: ";
            !read  *, n;                             !  User size input
            !print *, n                              !  Echo size
            
            !if ( n > MAX ) then                     !  Check validity
            !    stop "Vector size exceeds dimension, ERROR"; ! Abort run
            !end if
            
            !  Input vector contents
            !print *,"Enter the vector: "; 
            !read *, (x(i), i=1,n)              !  Read values
        
            vector_norm = sqrt ( sum ( x(:n)*x(:n) ))   ! L2 norm
            !L1_norm = sum ( abs(x(:n)) )/n          ! L1 norm 
        
            !write (6,'(" The L2 norm of [")',advance='no') ! Print header
        !   Loop over vector contents, (except last)
            !do i = 1, n-1
            !    write (6, '(f10.4, ",")', advance = 'no') x(i) ! List values
            !end do 
        
            !  Now last term
            !write (6, '(f10.4, " ] is: ")', advance = 'no') x(n) ! Last value
            
            !print *, L2_norm, "."  ;                    !  Print L2
            !print *, "The L1 norm is: ", L1_norm, "." ; !  Print L1
           
    end function vector_norm
    !********************************************************************************
    ! This subroutine is useful to expand integer arrays.
    !
    !  Array, x, will be allocated if the requested dimension is 1 (i.e., n=1).
    !  Array, x, will be expanded to the requested dimension, n, if (n > dim(x)).
    !
    ! For example, with the input
    !
    !  n = 7
    !  x = [9,4,2,1]
    !
    ! this subroutine will return
    !
    !  x = [9,4,2,1,0,0,0].
    !
    ! Note: If n=1, this subroutine takes is as an initialization, and returns
    !
    !  x = [0].
    !
    ! Note: So, this subroutine can only expand an interger array. It does not
    !       shrink an array. If n < size(x), then it is considered as an error and
    !       stop. If you want, I believe you can implement it.
    !
    !********************************************************************************
    subroutine my_alloc_int_ptr(x,n)

        implicit none
      
        integer, intent(in) :: n
      
        integer, dimension(:), pointer :: x
        integer, dimension(:), pointer :: temp
      
        integer :: i
      
        ! Error if n is negative.
      
        if (n <= 0) then
            write(*,*) "my_alloc_int_ptr received non-positive dimension. Stop."
            stop
        endif
      
        ! Shirinking an array is not implemented... Sorry.
        if ( n < size(x) ) then
            write(*,*) "my_alloc_int_ptr received a smaller dimension. Not implemented. Stop."
            stop
        endif
      
        ! If dimension 1, just allocate and return.
      
        if (n==1) then
            if (associated(x)) nullify(x)
            allocate(x(1))
            return
        endif
      
        ! If reallocation (i.e., n > size(x), create a pointer with a target with the requested dimension.
      
        allocate(temp(n))
        temp = 0
      
        ! Copy the existing data: e.g., for n=7 and x=[9,4,2,1] -> we construct temp=[9,4,2,1,0,0,0].
      
        do i = 1, size(x)
            temp(i) = x(i)
        end do
      
        ! Re-assign the pointer: x=[9,4,2,1,0,0,0].
            x => temp
      
        return
      
        end subroutine my_alloc_int_ptr
        !********************************************************************************
end module module_ccfv_gradient