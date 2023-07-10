module module_sparse_matrix
    ! This is a module for building and manipulating compressed sparse matrices specific to the linear solver 
    
    implicit none

    public :: build_A_BCSM ! builds A matrix in a Block Compressed Sparse Matrix (BCSM) format

contains
    
    subroutine build_A_BCSM(ncells,cell,jac,V,C,R,nnz)
        ! Takes in a cell structure and a jacobian structure and creates a corresponding A matrix using the Yale meethod:
        ! https://en.wikipedia.org/wiki/Sparse_matrix
        !
        ! Note : Requires Fortran Standard Library
        ! Installation instructions: https://github.com/fortran-lang/stdlib#getting-started

        use module_common_data, only : p2
        use module_ccfv_data_grid, only : cc_data_type, jacobian_type
        implicit none 

        integer, intent(in) :: ncells
        type(cc_data_type), dimension(ncells), intent(in) :: cell
        type(jacobian_type), dimension(ncells), intent(in) :: jac

        real(p2), dimension(:,:,:), allocatable, intent(out)  :: V   ! Values (5x5 block matrix) plus corresponding index
        integer, dimension(:),  allocatable, intent(out)       :: C   ! Column index of each value
        integer,dimension(ncells+1), intent(out) :: R   ! Start index of each new row
        integer, INTENT(OUT), optional :: nnz

        integer :: i, j, length

        R(1) = 1 ! Row 1 starts at 1
        do i = 2,ncells + 1
            R(i) = R(i-1) + 1 + cell(i-1)%nnghbrs ! Start of row(i) = row(i-1) start point + 1 (diagonal term) + # of neighbors
        end do
        ! R(ncells + 1) = R(ncells) + 1 ! "Dummy" variable marks end of last row (start of "ghost" next row)
        nnz = R(ncells+1) - 1

        allocate(V(5,5,nnz))
        allocate(C(    nnz))

        do i = 1,ncells
            ! sort the index of the cell neighbors and i and stores them in C:
            call insertion_sort_index( (/ cell(i)%nghbr, i /) , C(R(i) : (R(i+1)-1)) ) 
            length = R(i+1)-R(i)
            do j = R(i),(R(i+1)-1)
                if (length == C(j)) then
                    V(:,:,j) = jac(i)%diag(:,:)
                    C(j) = i
                else
                    V(:,:,j) = jac(i)%off_diag(:,:,C(j))
                    C(j) = cell(i)%nghbr(C(j))
                end if
            end do
        end do
    end subroutine build_A_BCSM

    subroutine A_times_P(ncells,ngroups,nnz,V,C,R,Restrict,productV,productC,productR)
        use module_common_data, only : p2, zero
        ! This subroutine performes the first half (AP) of RAP (where P = R^T).  In preparation for multiplication by R
        ! the resulting matrix is stored in Block Compressed Sparse COLUMN (BCSC) form.  This allows for full columns to be easily 
        ! accessed when it comes to multiply time R.  This is (currently) the only array that is stored this way.  

        implicit none

        integer, intent(in) :: ncells
        integer, intent(in) :: ngroups
        integer, intent(in) :: nnz
        real(p2), dimension(5,5,nnz), intent(in) :: V
        integer, dimension(nnz), intent(in) :: C
        integer, dimension(ncells+1), intent(in) :: R
        integer, dimension(ngroups,ncells), intent(in) :: Restrict

        real(p2), dimension(:,:,:), allocatable, intent(out) :: productV
        integer, dimension(:), allocatable, INTENT(OUT) :: productR
        integer, dimension(ngroups+1), INTENT(OUT) :: productC

        integer :: nnz_prime, i, j, k, os
        logical :: added
        real(p2), dimension(5,5) :: intermediateProduct
        character(len=80) :: filenameV, filenameC, filenameR
        filenameV = "A_V.dat"
        filenameC = "A_C.dat"
        filenameR = "A_R.dat"

        nnz_prime = 0
        ! count the number of nonzero entries in the product
        do k = 1,ngroups
            row_mid : do i = 1,ncells
                do j = R(i),(R(i+1)-1)
                    ! see if any of the values from column k of P match to any of the values of row i of A.  Since it's any not
                    ! all, this is essentially an OR gate.  Which means once we find one we can cycle the row.
                    if (Restrict(k,C(j)) == 1) then ! Restrict (j,i) = P(i,j)
                        nnz_prime = nnz_prime + 1
                        cycle row_mid
                    end if
                end do
            end do row_mid
        end do

        allocate(productV(5,5,nnz_prime))
        allocate(productR(    nnz_prime))
        
        ! **************************************************************
        ! BEGIN DEBUG CODE
        ! **************************************************************
        ! open(unit=8, file=filenameV, status="replace", iostat=os)   
        ! do i = 1,nnz
        !     do k = 1,5
        !         write(8,*) V(k,:,i)
        !     end do
        ! end do
        ! close(8)

        ! open(unit=8, file=filenameC, status="replace", iostat=os)   
        ! do i = 1,nnz
        !         write(8,*) C(i)
        ! end do
        ! close(8)


        ! open(unit=8, file=filenameR, status="replace", iostat=os)   
        ! do i = 1,ncells+1
        !         write(8,*) R(i)
        ! end do
        ! close(8)

        ! write(*,*) "NNZ = ", nnz
        ! write(*,*) "NNZ_prime = ", nnz_prime


        
        ! **************************************************************
        ! END DEBUG CODE
        ! **************************************************************

        nnz_prime = 1 ! re-init
        productC(1) = nnz_prime
        productV = zero

        do k = 1,ngroups
            do i = 1,ncells
                ! if (nnz_prime == 916) then
                !     write(*,*) productV(1,1,nnz_prime)
                ! end if
                added = .false.
                do j = R(i),(R(i+1)-1)
                    intermediateProduct = V(:,:,j) * Restrict(k,C(j))
                    if ( any(intermediateProduct(:,:) /= zero) ) then
                        productV(:,:,nnz_prime) = productV(:,:,nnz_prime) + intermediateProduct
                        productR(nnz_prime) = i
                        added = .true.
                    end if
                end do
                ! if (any(productV(:,:,nnz_prime) /= zero)) then
                !     write(*,*) productV(1,1,nnz_prime)
                ! end if
                if (added) nnz_prime = nnz_prime + 1
            end do
            productC(k+1) = nnz_prime
        end do
        ! Note: nnz_prime is now one higher than it should be.  That is fine since it isn't being passed up.  If that changes we 
        ! will need to add:
        ! nnz_prime = nnz_prime - 1

    end subroutine A_times_P


    subroutine R_A_P(ncells,ngroups,nnz,Restrict,V,C,R,RAP_V,RAP_C,RAP_R,nnz_prime_final)
        ! This subroutine computes the restriction matrix A' = RAP for algebraic multi grid and stores it using BCSM.
        ! Note: since piecewise injection is being used P = (R^T), so only R is needed (indices are just swapped).
        use module_common_data , only : p2, zero

        implicit none
        integer, intent(in) :: ncells                                   ! # of cells for the coarse mesh
        integer, intent(in) :: ngroups                                  ! # of groups on the restricted level
        integer, intent(in) :: nnz                                      ! # of nonzero entries in the fine A matrix
        integer, dimension(ngroups,ncells), intent(in) :: Restrict      ! Restriction matrix
        real(p2), dimension(5,5,nnz), intent(in) :: V                   ! Block Sparse Compressed Matrix (BSCM) values of A matrix
        integer, dimension(nnz), intent(in) :: C                        ! BCSM Columns of A matrix
        integer, dimension(ncells+1), intent(in) :: R                   ! BCSM Rows of A matrix

        real(p2), dimension(:,:,:), allocatable, INTENT(OUT) :: RAP_V   ! BCSM Values of restricted A matrix
        integer, dimension(:), allocatable, INTENT(OUT) :: RAP_C        ! BCSM Columns of restricted A matrix
        integer, dimension(ngroups + 1), INTENT(OUT) :: RAP_R           ! BCSM Rows of restricted A matrix
        integer, optional, intent(out) :: nnz_prime_final               ! (Optional) # of nonzero values in restricted A matrix

        real(p2), dimension(:,:,:), allocatable :: AP_V                 ! BCSM Values of intermediate A*(R^T) matrix
        integer, dimension(:), allocatable :: AP_R                      ! BCSM Columns of intermediate A*(R^T) matrix
        integer, dimension(ngroups + 1):: AP_C                          ! BCSM Rows of intermediate A*(R^T) matrix
        
        integer :: i, j, k, nnz_prime
        
        nnz_prime = 0
        call A_times_P(ncells,ngroups,nnz,V,C,R,Restrict,AP_V,AP_C,AP_R)

        do k = 1,ngroups
            mid_column : do i = 1,ngroups
                do j = AP_C(i),(AP_C(i+1)-1)
                    if (Restrict(k,AP_R(j)) == 1) then
                        nnz_prime = nnz_prime + 1
                        cycle mid_column
                    end if
                end do
            end do mid_column
        end do

        allocate(RAP_V(5,5,nnz_prime))
        allocate(RAP_C(    nnz_prime))
        RAP_V = zero
        ! reset nnz_prime (to 1)
        RAP_R(1) = 1
        nnz_prime = 1
        do k = 1,ngroups
            do i = 1,ngroups
                do j = AP_C(i),(AP_C(i+1)-1)
                    RAP_V(:,:,nnz_prime) = RAP_V(:,:,nnz_prime) + AP_V(:,:,j) * Restrict(k,AP_R(j))
                end do 
                RAP_C(nnz_prime) = i
                if (any(RAP_V(:,:,nnz_prime) /= zero)) nnz_prime = nnz_prime + 1
            end do
            RAP_R(k+1) = nnz_prime
        end do

        deallocate(AP_V)
        deallocate(AP_R)

        if (present(nnz_prime_final)) nnz_prime_final = nnz_prime


    end subroutine R_A_P 

    subroutine insertion_sort_index(x,index)
        ! This routine uses insertion sort sort the input vector x and returns the output vector index
        ! which is the index needed to sort x.
        ! Insertion sort is O(n^2) which is generally inefficient.  However its simplicity allows it 
        ! to be faster than most O(n log n) methods for small arrays.  Since the incoming array will 
        ! just be the cell neighbors for the coarse mesh n <= 6 unless I get around to adding support
        ! for general polyhedral cells (please no...)
        !
        ! ref: https://en.wikipedia.org/wiki/Insertion_sort

        integer, dimension(:), intent(in)  :: x
        integer, dimension(:), intent(out) :: index
        
        integer, dimension(:), allocatable :: x_internal
        integer i, j, length, x_j, index_j
        length = size(x)
        allocate(x_internal(length))
        x_internal = x

        index = (/(i,i=1,length)/)
        
        do i = 1,length
            j = i
            inner : do 
                if ( (j <= 1) ) then
                    exit inner
                else if ( x_internal(j-1) < x_internal(j)) then ! it doesn't like evaluating this line with .or. 
                    !if j = 1 (array bounds)
                    ! I know this is stupid, but I don't want to debug it right now...
                    exit inner
                end if
                x_j = x_internal(j)
                x_internal(j) = x_internal(j-1)
                x_internal(j-1) = x_j
                index_j = index(j)
                index(j) = index(j-1)
                index(j-1) = index_j
                j = j-1
            end do inner
        end do
        deallocate(x_internal)
    end subroutine insertion_sort_index

    subroutine sparseblock_times_vectorblock(ncells,V,C,R,x,b)
        use module_common_data , only : p2, zero

        implicit none
        
        integer, intent(in)                    :: ncells
        real(p2), dimension(:,:,:), intent(in) :: V
        integer, dimension(:), intent(in)      :: C
        integer, dimension(ncells+1), intent(in)      :: R
        real(p2), dimension(5,ncells), intent(in)     :: x

        real(p2), dimension(5,ncells), intent(out)    :: b

        real(p2), dimension(5,5)  :: currentV
        real(p2), dimension(5)  :: currentX

        integer :: i,j,ii
        ! real(p2), dimension(5) :: product_test
        b = zero
        do i = 1,ncells
            do j = R(i),(R(i+1)-1)
                ! currentV(:,:) = V(:,:,j)
                ! currentX(:)   = X(:,C(j))
                ! write(*,*) "V(:,:,", j ,"):"
                ! do ii = 1,5
                !     write(*,*) V(ii,:,j)
                ! end do
                ! write(*,*) "X(:,C(",j,"))"
                ! write(*,*) x(:,C(j))
                b(:,i) = b(:,i) + matmul(V(:,:,j),x(:,C(j)));
                ! write(*,*) matmul(V(:,:,j),x(:,C(j)))
            end do
        end do
        
    end subroutine sparseblock_times_vectorblock

    subroutine vector_times_sparse(ncells,V,C,R,x,b)
        use module_common_data , only : p2

        implicit none
        
        integer, intent(in)                    :: ncells
        real(p2), dimension(:), intent(in) :: V
        integer, dimension(:), intent(in)      :: C
        integer, dimension(ncells+1), intent(in)      :: R
        real(p2), dimension(ncells), intent(in)     :: x

        real(p2), dimension(ncells), intent(out)    :: b

        integer :: i, j

        do i = 1,ncells
            do j = R(i),(R(i+1)-1)
                b(C(j)) = b(C(j)) +  x(i) * V(j)
            end do
        end do

    end subroutine vector_times_sparse
    subroutine destroy_A(V,C)
        use module_common_data , only : p2
        ! Deallocates any allocated value and column arrays
        real(p2), dimension(:,:,:),optional, allocatable, intent(out)  :: V   ! Values (5x5 block matrix) plus corresponding index
        integer, dimension(:),optional, allocatable, intent(out)       :: C   ! Column index of each value

        if (allocated(V)) deallocate(V)
        if (allocated(C)) deallocate(C)

    end subroutine
end module module_sparse_matrix