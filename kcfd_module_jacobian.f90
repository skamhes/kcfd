module module_jacobian
    implicit none
    
    public :: compute_jacobian ! compute the jacobian
    public :: gewp_solve       ! Gauss elimination for inverting diagonal blocks
contains
    
    subroutine compute_jacobian
        
        use module_common_data          , only : p2, zero, nb, bc_type
        use module_flux_jac_interface   , only : interface_jac
        use module_ccfv_data_grid       , only : nfaces, face, ncells, cell, face_nrml, face_nrml_mag, jac, &
                                                 kth_nghbr_of_1, kth_nghbr_of_2, bound
        use module_ccfv_data_soln       , only : w, u, u2w, dtau
        use module_bc_states            , only : get_right_state
        
        implicit none
        ! Local Vars
        integer                     :: c1, c2, i, k, ib, idestat, j
        real(p2), dimension(3)      :: unit_face_nrml, bface_centroid
        real(p2), dimension(5)      :: u1, ub, wb
        real(p2), dimension(5,5)    :: dFnduL, dFnduR
        real(p2)                    :: face_mag
        ! Initialize jacobian terms
        do i = 1,ncells
            jac(i)%diag = zero
            jac(i)%off_diag = zero
            jac(i)%diag_inv = zero
        end do
        ! Loop Faces
        do i = 1,nfaces
            c1 = face(1,i)
            c2 = face(2,i)

            unit_face_nrml = face_nrml(1:3,i)
            face_mag       = face_nrml_mag(i)
            ! Compute the flux Jacobian for given w1 and w2
            call interface_jac( w(:,c1), w(:,c2), unit_face_nrml, dFnduL, dFnduR)

            ! Add to diagonal term of C1
            jac(c1)%diag            = jac(c1)%diag            + dFnduL * face_mag
            ! get neighbor index k for cell c1
            k = kth_nghbr_of_1(i)
            ! add to off diagonal neighbor k for cell c1
            jac(c1)%off_diag(:,:,k) = jac(c1)%off_diag(:,:,k) + dFnduR * face_mag

            ! Subtract terms from c2
            jac(c2)%diag            = jac(c2)%diag            - dFnduR * face_mag
            k = kth_nghbr_of_2(i)
            jac(c2)%off_diag(:,:,k) = jac(c2)%off_diag(:,:,k) - dFnduL * face_mag
            
        end do
        
        ! Now we gotta do the boundary faces
        bound_loop : do ib = 1,nb
            bfaces_loop : do i = 1,bound(ib)%nbfaces
                c1 = bound(ib)%bcell(i)
                u1 = u(1:5,c1)
                bface_centroid = bound(ib)%bface_center(:,i)
                unit_face_nrml = bound(ib)%bface_nrml(:,i)
                face_mag       = bound(ib)%bface_nrml_mag(i)
                call get_right_state(bface_centroid(1), bface_centroid(2), bface_centroid(3), &
                                     u1, unit_face_nrml, bc_type(ib), ub)
                wb = u2w(ub)
                call interface_jac( w(:,c1), wb, unit_face_nrml, dFnduL, dFnduR)
                ! We only have a diagonal term to add
                jac(c1)%diag            = jac(c1)%diag            + dFnduL * face_mag
            end do bfaces_loop
        end do bound_loop

        !--------------------------------------------------------------------------------
        ! Add pseudo time term vol/dtau to the diagonals of the diagonal blocks
        ! to form: Jacobian = V/dtau + dRes1/dU, where V/dtau is a global diagonal
        ! matrix having vol(i)/dtau(i) for each node, and Res1 is the first-order
        ! residual.
        do i = 1,ncells
            do k = 1,5
                jac(i)%diag(k,k) = jac(i)%diag(k,k) + cell(i)%vol/dtau(i)
            end do
        end do

        ! Invert the diagonal block matrices using Gauss elimination with pivoting
        !  Note: gewp_solve() actually solves a linear system, Ax=b, but here
        !        we use it only to obtain the inverse of A. So, b is a dummy.
        cell_inverse_loop : do i = 1,ncells
            !jac(i)%diag_inv = zero
            idestat = 0
            !                A                 dim  A^{-1}           error check
            call gewp_solve( jac(i)%diag(:,:), 5  , jac(i)%diag_inv, idestat    )
             !  Report errors
            if (idestat/=0) then
                write(*,*) " Error in inverting the diagonal block... Stop"
                write(*,*) "  Cell number = ", i
                do k = 1, 5
                    write(*,'(12(es8.1))') ( jac(i)%diag(k,j), j=1,5 )
                end do
                stop
            endif
        end do cell_inverse_loop
        
    end subroutine compute_jacobian



    !****************************************************************************
    !* ------------------ GAUSS ELIMINATION WITH PIVOTING ---------------------
    !*
    !*  This computes the inverse of an (nm)x(nm) matrix "ai" and also
    !*  computes the solution to a given lienar system.
    !*
    !*  IN :       ai = An (nm)x(nm) matrix whoise inverse is sought.
    !*             bi = A vector of (nm): Right hand side of the linear sytem
    !*             nm = The size of the matrix "ai"
    !*
    !* OUT :
    !*            sol = Solution to the linear system: ai*sol=bi
    !*        inverse = the inverse of "ai".
    !*       idetstat = 0 -> inverse successfully computed
    !*                  1 -> THE INVERSE DOES NOT EXIST (det=0).
    !*                  2 -> No unique solutions exist.
    !*****************************************************************************
    subroutine gewp_solve(ai,nm, inverse,idetstat)

        implicit none
      
        integer , parameter ::    p2 = selected_real_kind(15) ! Double precision
        real(p2), parameter ::  zero = 0.0_p2
        real(p2), parameter ::   one = 1.0_p2
      
        integer ,                   intent( in) :: nm
        real(p2), dimension(nm,nm), intent( in) :: ai
      
        real(p2), dimension(nm,nm), intent(out) :: inverse
        integer ,                   intent(out) :: idetstat
      
        real(p2), dimension(nm,nm+1) :: a
        real(p2), dimension(nm)      :: x
        integer , dimension(nm)      :: nrow
        integer                      :: I,J,K,pp,m
      
        do m = 1, nm
            !*****************************************************************************
            !* Set up the matrix a
            !*****************************************************************************
            
            do J=1,nm
                do I=1,nm
                a(I,J) = ai(I,J)
                end do
            end do
        
            do k=1,nm
                a(k,nm+1)=zero; nrow(k)=k
            end do
            a(m,nm+1)=one
        
            !*****************************************************************************
            !* HONA IKOKA..... 
            !*****************************************************************************
            do j=1,nm-1
            !*****************************************************************************
            !* FIND SMALLEST pp FOR a(pp,j) IS MAXIMUM IN JTH COLUMN.
            !***************************************************************************** 
                call findmax(nm,j,pp,a,nrow)
                !*****************************************************************************
                !* IF a(nrow(p),j) IS zero, THERE'S NO UNIQUE SOLUTIONS      
                !*****************************************************************************
                if (abs(a(nrow(pp),j)) < epsilon(one)) then
                    write(6,*) 'THE INVERSE DOES NOT EXIST.'
                    idetstat = 1
                    return
                endif
                !*****************************************************************************
                !* IF THE MAX IS NOT A DIAGONAL ELEMENT, SWITCH THOSE ROWS       
                !*****************************************************************************
                if (nrow(pp) .ne. nrow(j)) then
                    call switch(nm,j,pp,nrow)
                else
                endif  
                !*****************************************************************************
                !* ELIMINATE ALL THE ENTRIES BELOW THE DIAGONAL ONE
                !***************************************************************************** 
                call eliminate_below(nm,j,a,nrow)
        
            end do
            !*****************************************************************************
            !* CHECK IF a(nrow(N),N)=0.0 .
            !*****************************************************************************
            if (abs(a(nrow(nm),nm)) < epsilon(one)) then
                write(6,*) 'NO UNIQUE SOLUTION EXISTS!'
                idetstat = 2
                return
            else
            endif
            !*****************************************************************************
            !* BACKSUBSTITUTION!
            !*****************************************************************************
            call backsub(nm,x,a,nrow)
            !*****************************************************************************
            !* STORE THE SOLUTIONS, YOU KNOW THEY ARE INVERSE(i,m) i=1...
            !*****************************************************************************
            do i=1,nm
                inverse(i,m)=x(i)
            end do
            !*****************************************************************************
        end do
      
        idetstat = 0
      
        return
      
        !*****************************************************************************
    end subroutine gewp_solve
      
      !*****************************************************************************
      !* Four subroutines below are used in gewp_solve() above.
      !*****************************************************************************
      !* FIND MAXIMUM ELEMENT IN jth COLUMN 
      !***************************************************************************** 
            subroutine findmax(nm,j,pp,a,nrow)
      
            implicit none
      
            integer , parameter   :: p2 = selected_real_kind(15) ! Double precision
            integer , intent( in) :: nm
            real(p2), intent( in) :: a(nm,nm+1)
            integer , intent( in) :: j,nrow(nm)
            integer , intent(out) :: pp
            real(p2)              :: max
            integer               :: i
      
                  max=abs(a(nrow(j),j)); pp=j
      
                 do i=j+1,nm
      
                   if (max < abs(a(nrow(i),j))) then
      
                        pp=i; max=abs(a(nrow(i),j))
      
                   endif
      
                 end do
      
            return
      
            end subroutine findmax
      !*****************************************************************************
      !* SWITCH THOSE ROWS       
      !*****************************************************************************
            subroutine switch(nm,j,pp,nrow)
      
            implicit none
      
            integer, intent(   in) :: nm,j,pp
            integer, intent(inout) :: nrow(nm)
            integer                :: ncopy
      
            if (nrow(pp).ne.nrow(j)) then
      
               ncopy=nrow(j)
               nrow(j)=nrow(pp)
               nrow(pp)=ncopy
      
            endif
      
            return
      
            end subroutine switch
      !*****************************************************************************
      !* ELIMINATE ALL THE ENTRIES BELOW THE DIAGONAL ONE
      !*(Give me j, the column you are working on now)
      !***************************************************************************** 
            subroutine eliminate_below(nm,j,a,nrow)
      
            implicit none
      
            integer , parameter     :: p2 = selected_real_kind(15) ! Double precision
            real(p2), parameter     :: zero = 0.0_p2
            integer , intent(   in) :: nm
            real(p2), intent(inout) :: a(nm,nm+1)
            integer , intent(   in) :: j,nrow(nm)
            real(p2)                :: m
            integer                 :: k,i
      
            do i=j+1,nm
      
              m=a(nrow(i),j)/a(nrow(j),j)
              a(nrow(i),j)=zero
      
                do k=j+1,nm+1
                  a(nrow(i),k)=a(nrow(i),k)-m*a(nrow(j),k)
                end do
      
            end do
      
            return
      
            end subroutine eliminate_below
      !*****************************************************************************
      !* BACKSUBSTITUTION!
      !*****************************************************************************
            subroutine backsub(nm,x,a,nrow)
      
            implicit none
      
            integer , parameter   :: p2 = selected_real_kind(15) ! Double precision
            real(p2), parameter   :: zero = 0.0_p2
      
            integer , intent( in) :: nm
            real(p2), intent( in) :: a(nm,nm+1)
            integer , intent( in) :: nrow(nm)
            real(p2), intent(out) :: x(nm)
            real(p2)              :: sum
            integer               :: i,k
      
            x(nm)=a(nrow(nm),nm+1)/a(nrow(nm),nm)
      
            do i=nm-1,1,-1
      
               sum=zero
      
                 do k=i+1,nm
      
                    sum=sum+a(nrow(i),k)*x(k)
      
                 end do
      
            x(i)=(a(nrow(i),nm+1)-sum)/a(nrow(i),i)
      
            end do
      
            return
      
            end subroutine backsub
      !*********************************************************************
      
end module module_jacobian