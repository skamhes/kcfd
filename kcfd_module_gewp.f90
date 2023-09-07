module module_gewp
    implicit none
    
    public :: gewp_solve       ! Gauss elimination for inverting diagonal blocks


contains

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

end module module_gewp