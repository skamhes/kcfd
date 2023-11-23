!********************************************************************************
! Educationally-Designed Unstructured 3D (EDU3D) Code
!
!  --- EDU3D Euler
!
!
! This module containes subroutines that computes numerical fluxes for the
! 3D Euler equations. The following flux functions are implemented: 
!
!  - Roe
!  - HLL
!  - Rusanov
!  - RHLL
!
! Note: see edu3d_module_ddt5.f90 for details on how automatic differentiation
!       works.
!
!
!
!        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!
! the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!
! 
! This F90 program is written and made available for an educational purpose.
!
! This file may be updated in future.
!
! Katate Masatsuka, October 2017. http://www.cfdbooks.com
!********************************************************************************

 module flux_functions_ddt

 !This module contains the following subroutines:

  public ::            roe_ddt
  public ::        rusanov_ddt
  public ::            hll_ddt
  public ::           rhll_ddt
  public ::  viscous_alpha_ddt

 contains

!********************************************************************************
!* -- 3D Roe's Flux Function and Jacobian --
!*
!* NOTE: This version does not use any tangent vector.
!*       See "I do like CFD, VOL.1" about how tangent vectors are eliminated.
!*
!* This subroutine computes the Roe flux for the Euler equations
!* in the direction, njk=[nx,ny,nz].
!*
!* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
!* Schemes, Journal of Computational Physics, 43, pp. 357-372.
!*
!* Conservative form of the Euler equations:
!*
!*     dU/dt + dF/dx + dG/dy + dH/dz = 0
!*
!* This subroutine computes the numerical flux for the flux in the direction,
!* njk=[nx,ny,nz]:
!*
!*     Fn = F*nx + G*ny + H*nz = | rho*qn          |
!*                               | rho*qn*u + p*nx |
!*                               | rho*qn*v + p*ny |
!*                               | rho*qn*w + p*nz |
!*                               | rho*qn*H        |    (qn = u*nx + v*ny + w*nz)
!*
!* The Roe flux is implemented in the following form:
!*
!*   Numerical flux = 1/2 [ Fn(UR) + Fn(UL) - |An|dU ], 
!*
!*  where
!*
!*    An = dFn/dU,  |An| = R|Lambda|L, dU = UR - UL.
!*
!* The dissipation term, |An|dU, is actually computed as
!*
!*     sum_{k=1,4} |lambda_k| * (LdU)_k * r_k,
!*
!* where lambda_k is the k-th eigenvalue, (LdU)_k is the k-th wave strength,
!* and r_k is the k-th right-eigenvector evaluated at the Roe-average state.
!*
!* Note: The 4th component is a combined contribution from two shear waves.
!*       They are combined to eliminate the tangent vectors.
!*       So, (LdU)_4 is not really a wave strength, and
!*       r_4 is not really an eigenvector.
!*       See "I do like CFD, VOL.1" about how tangent vectors are eliminated.
!*
!* Note: In the code, the vector of conserative variables are denoted by uc.
!*
!* ------------------------------------------------------------------------------
!*  Input: ucL(1:5) =  Left state (rhoL, rhoL*uL, rhoL*vL, rhoL*wR, rhoL*EL)
!*         ucR(1:5) = Right state (rhoR, rhoL*uR, rhoL*vR, rhoL*wR, rhoL*ER)
!*         njk(1:3) = unit face normal vector (nx, ny, nz), pointing from Left to Right.
!*
!*           njk
!*  Face normal ^   o Right data point
!*              |  .
!*              | .
!*              |. 
!*       -------x-------- Face
!*             .                 Left and right states are
!*            .                   1. Values at data points for 1st-order accuracy
!*           .                    2. Extrapolated values at the face midpoint 'x'
!*          o Left data point        for 2nd/higher-order accuracy.
!*
!*
!* Output:  num_flux(1:5) = the numerical flux vector
!*          dFdU(1:5,1:5) = flux Jacobian matrix
!*                    wsn = maximum wave speed (eigenvalue)
!*
!* ------------------------------------------------------------------------------
!*
!* Note: This function is a ddt-version, which means that each variable carries
!*       its derivatives, and that the resulting flux "numerical_flux" will have
!*       its derivatives in "numerical_flux%df".
!*
!* Note: This subroutine has been prepared for an educational purpose.
!*       It is not at all efficient. Think about how you can optimize it.
!*       One way to make it efficient is to reduce the number of local variables,
!*       by re-using temporary variables as many times as possible.
!*
!* Note: Please let me know if you find bugs. I'll greatly appreciate it and
!*       fix the bugs.
!*
!* Katate Masatsuka, November 2012. http://www.cfdbooks.com
!********************************************************************************
 subroutine roe_ddt(ucL, ucR, njk, num_flux,dFdU,wsn)

  use derivative_data_df5
 
  implicit none
  integer , parameter :: p2 = selected_real_kind(15) ! Double precision
 
 !Input
  type(derivative_data_type_df5), dimension(5), intent( in) :: ucL
  type(derivative_data_type_df5), dimension(5), intent( in) :: ucR
  real(p2)                      , dimension(3), intent( in) :: njk
 
 !Output
  real(p2), dimension(5)  , intent(out) :: num_flux !Numerical viscous flux
  real(p2), dimension(5,5), intent(out) :: dFdU     !Numerical viscous flux Jacobian
  real(p2),                 intent(out) :: wsn      ! Max wave speed
 
 !Some constants
  real(p2) ::    zero = 0.00_p2
  real(p2) ::     one = 1.00_p2
  real(p2) ::     two = 2.00_p2
  real(p2) ::    half = 0.50_p2
  real(p2) ::   gamma = 1.40_p2                 ! Ratio of specific heats
 
 !Local variables
 !
 !            L = Left
 !            R = Right
 ! No subscript = Roe average
 
  type(derivative_data_type_df5), dimension(5) :: numerical_flux !Numerical flux in ddt
  real(p2)                       :: nx, ny, nz             ! Normal vector components
 
  type(derivative_data_type_df5) :: uL, uR, vL, vR, wL, wR ! Velocity components.
  type(derivative_data_type_df5) :: rhoL, rhoR, pL, pR     ! Primitive variables.
  type(derivative_data_type_df5) :: qnL, qnR               ! Normal velocities
  type(derivative_data_type_df5) :: aL, aR, HL, HR         ! Speed of sound, Total enthalpy
  type(derivative_data_type_df5), dimension(5)   :: fL     ! Physical flux evaluated at ucL
  type(derivative_data_type_df5), dimension(5)   :: fR     ! Physical flux evaluated at ucR
 
  type(derivative_data_type_df5) :: RT                     ! RT = sqrt(rhoR/rhoL)
  type(derivative_data_type_df5) :: rho,u,v,w,H,a,qn       ! Roe-averages
 
  type(derivative_data_type_df5) :: drho, dqn, dp          ! Differences in rho, qn, p, e.g., dp=pR-pL
  type(derivative_data_type_df5), dimension(4) :: LdU      ! Wave strengths = L*(UR-UL)
  type(derivative_data_type_df5) :: du, dv, dw             ! Velocity differences
  type(derivative_data_type_df5), dimension(4)   :: ws     ! Wave speeds
  type(derivative_data_type_df5), dimension(4)   :: dws    ! Width of a parabolic fit for entropy fix
  type(derivative_data_type_df5), dimension(5,4) :: R      ! Right-eigenvector matrix
  type(derivative_data_type_df5), dimension(5)   :: diss   ! Dissipation term
 
  type(derivative_data_type_df5) :: temp
  integer                        :: i, j
 
 ! Face normal vector (unit vector)
 
   nx = njk(1)
   ny = njk(2)
   nz = njk(3)
 
 !Primitive and other variables.
 
 !  Left state
 
     rhoL = ucL(1)
       uL = ucL(2)/ucL(1)
       vL = ucL(3)/ucL(1)
       wL = ucL(4)/ucL(1)
      qnL = uL*nx + vL*ny + wL*nz
       pL = (gamma-one)*( ucL(5) - half*rhoL*(uL*uL+vL*vL+wL*wL) )
       aL = ddt_sqrt(gamma*pL/rhoL)
       HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL+wL*wL)
 
 !  Right state
 
     rhoR = ucR(1)
       uR = ucR(2)/ucR(1)
       vR = ucR(3)/ucR(1)
       wR = ucR(4)/ucR(1)
      qnR = uR*nx + vR*ny + wR*nz
       pR = (gamma-one)*( ucR(5) - half*rhoR*(uR*uR+vR*vR+wR*wR) )
       aR = ddt_sqrt(gamma*pR/rhoR)
       HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR+wR*wR)
 
 !Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)
 
   fL(1) = rhoL*qnL
   fL(2) = rhoL*qnL * uL + pL*nx
   fL(3) = rhoL*qnL * vL + pL*ny
   fL(4) = rhoL*qnL * wL + pL*nz
   fL(5) = rhoL*qnL * HL
 
   fR(1) = rhoR*qnR
   fR(2) = rhoR*qnR * uR + pR*nx
   fR(3) = rhoR*qnR * vR + pR*ny
   fR(4) = rhoR*qnR * wR + pR*nz
   fR(5) = rhoR*qnR * HR
 
 !First compute the Roe-averaged quantities
 
 !  NOTE: See http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
 !        the Roe-averaged density.
 
     RT = ddt_sqrt(rhoR/rhoL)
    rho = RT*rhoL                                            !Roe-averaged density
      u = (uL + RT*uR)/(one + RT)                            !Roe-averaged x-velocity
      v = (vL + RT*vR)/(one + RT)                            !Roe-averaged y-velocity
      w = (wL + RT*wR)/(one + RT)                            !Roe-averaged z-velocity
      H = (HL + RT*HR)/(one + RT)                            !Roe-averaged total enthalpy
      a = ddt_sqrt( (gamma-one)*(H-half*(u*u + v*v + w*w)) ) !Roe-averaged speed of sound
     qn = u*nx + v*ny + w*nz                                 !Roe-averaged face-normal velocity
 
 !Wave Strengths
 
    drho = rhoR - rhoL !Density difference
      dp =   pR - pL   !Pressure difference
     dqn =  qnR - qnL  !Normal velocity difference
 
   LdU(1) = (dp - rho*a*dqn )/(two*a*a) !Left-moving acoustic wave strength
   LdU(2) =  drho - dp/(a*a)            !Entropy wave strength
   LdU(3) = (dp + rho*a*dqn )/(two*a*a) !Right-moving acoustic wave strength
   LdU(4) = rho                         !Shear wave strength (not really, just a factor)
 
 !Absolute values of the wave Speeds
 
   ws(1) = ddt_abs(qn-a) !Left-moving acoustic wave
   ws(2) = ddt_abs(qn)   !Entropy wave
   ws(3) = ddt_abs(qn+a) !Right-moving acoustic wave
   ws(4) = ddt_abs(qn)   !Shear waves
 
 ! Harten's Entropy Fix JCP(1983), 49, pp357-393. This is typically applied
 ! only for the nonlinear fields (k=1 and 3), but here it is applied to all
 ! for robustness, avoiding vanishing wave speeds by making a parabolic fit
 ! near ws = 0 for all waves.
 
   do i = 1, 4
    dws(i) = 0.2_p2*a
    if ( ws(i) < dws(i) ) ws(i) = half * ( ws(i)*ws(i)/dws(i)+dws(i) )
   end do
 
 !Right Eigenvectors
 !Note: Two shear wave components are combined into one, so that tangent vectors
 !      are not required. And that's why there are only 4 vectors here.
 !      See "I do like CFD, VOL.1" about how tangent vectors are eliminated.
 
 ! Left-moving acoustic wave
   R(1,1) = one    
   R(2,1) = u - a*nx
   R(3,1) = v - a*ny
   R(4,1) = w - a*nz
   R(5,1) = H - a*qn
 
 ! Entropy wave
   R(1,2) = one
   R(2,2) = u
   R(3,2) = v 
   R(4,2) = w
   R(5,2) = half*(u*u + v*v + w*w)
 
 ! Right-moving acoustic wave
   R(1,3) = one
   R(2,3) = u + a*nx
   R(3,3) = v + a*ny
   R(4,3) = w + a*nz
   R(5,3) = H + a*qn
 
 ! Two shear wave components combined into one (wave strength incorporated).
   du = uR - uL
   dv = vR - vL
   dw = wR - wL
   R(1,4) = zero
   R(2,4) = du - dqn*nx
   R(3,4) = dv - dqn*ny
   R(4,4) = dw - dqn*nz
   R(5,4) = u*du + v*dv + w*dw - qn*dqn
 
 !Dissipation Term: |An|(UR-UL) = R|Lambda|L*dU = sum_k of [ ws(k) * R(:,k) * L*dU(k) ]
 
  diss(:) = ws(1)*LdU(1)*R(:,1) + ws(2)*LdU(2)*R(:,2) &
          + ws(3)*LdU(3)*R(:,3) + ws(4)*LdU(4)*R(:,4)
 
 ! This is the numerical flux: Roe flux = 1/2 *[  Fn(UL)+Fn(UR) - |An|(UR-UL) ]
 
   numerical_flux = half * (fL + fR - diss)
 
 !--------------
 ! Output
 
    ! Normal max wave speed
     temp = ddt_abs(qn) + a
      wsn = temp%f
 
   do i = 1, 5
 
    !Numerical flux
     num_flux(i) = numerical_flux(i)%f
 
    do j = 1, 5
 
      !Flux derivative
       dFdU(i,j) = numerical_flux(i)%df(j)
 
    end do
 
   end do
 
  end subroutine roe_ddt
 !--------------------------------------------------------------------------------

  subroutine roe_low_mach_ddt(ucL, ucR, uR2L, uR2R, njk, num_flux,dFdU,wsn)

    use derivative_data_df5
   
    implicit none
    integer , parameter :: p2 = selected_real_kind(15) ! Double precision
   
   !Input
    type(derivative_data_type_df5), dimension(5), intent( in) :: ucL
    type(derivative_data_type_df5), dimension(5), intent( in) :: ucR
    real(p2)                                    , intent( in) :: uR2L, uR2R
    real(p2)                      , dimension(3), intent( in) :: njk
   
   !Output
    real(p2), dimension(5)  , intent(out) :: num_flux !Numerical viscous flux
    real(p2), dimension(5,5), intent(out) :: dFdU     !Numerical viscous flux Jacobian
    real(p2),                 intent(out) :: wsn      ! Max wave speed
   
   !Some constants
    real(p2) ::    zero = 0.00_p2
    real(p2) ::     one = 1.00_p2
    real(p2) ::     two = 2.00_p2
    real(p2) ::    half = 0.50_p2
    real(p2) ::   gamma = 1.40_p2                 ! Ratio of specific heats
   
   !Local variables
   !
   !            L = Left
   !            R = Right
   ! No subscript = Roe average
   
    type(derivative_data_type_df5), dimension(5) :: numerical_flux !Numerical flux in ddt
    real(p2)                       :: nx, ny, nz             ! Normal vector components
   
    type(derivative_data_type_df5) :: uL, uR, vL, vR, wL, wR ! Velocity components.
    type(derivative_data_type_df5) :: rhoL, rhoR, pL, pR     ! Primitive variables.
    type(derivative_data_type_df5) :: qnL, qnR               ! Normal velocities
    type(derivative_data_type_df5) :: aL, aR, HL, HR         ! Speed of sound, Total enthalpy
    type(derivative_data_type_df5), dimension(5)   :: fL     ! Physical flux evaluated at ucL
    type(derivative_data_type_df5), dimension(5)   :: fR     ! Physical flux evaluated at ucR
   
    type(derivative_data_type_df5) :: RT                     ! RT = sqrt(rhoR/rhoL)
    type(derivative_data_type_df5) :: rho,u,v,w,H,a,qn,p     ! Roe-averages
   
    type(derivative_data_type_df5) :: drho, dqn, dp          ! Differences in rho, qn, p, e.g., dp=pR-pL
    type(derivative_data_type_df5) :: drhou, drhov, drhow, drhoE ! Differences in conserved vars
    type(derivative_data_type_df5) :: absU
    type(derivative_data_type_df5) :: uprime, cprime, alpha, beta
    type(derivative_data_type_df5) :: uR2L_ddt, uR2R_ddt, uR2
    ! real(p2)                       :: uR2
    type(derivative_data_type_df5) :: cstar, Mstar, delu, delp
    type(derivative_data_type_df5), dimension(4) :: LdU      ! Wave strengths = L*(UR-UL)
    type(derivative_data_type_df5) :: du, dv, dw             ! Velocity differences
    type(derivative_data_type_df5), dimension(4)   :: ws     ! Wave speeds
    type(derivative_data_type_df5), dimension(4)   :: dws    ! Width of a parabolic fit for entropy fix
    type(derivative_data_type_df5), dimension(5,4) :: R      ! Right-eigenvector matrix
    type(derivative_data_type_df5), dimension(5)   :: diss   ! Dissipation term
   
    type(derivative_data_type_df5) :: temp
    
    type(derivative_data_type_df5), dimension(3,5) :: diss_cartesian
    real(p2) :: epsilon = 10e-05_p2
    integer                        :: i, j
   
   ! Face normal vector (unit vector)
   
     nx = njk(1)
     ny = njk(2)
     nz = njk(3)
   
   !Primitive and other variables.
   
   !  Left state
   
       rhoL = ucL(1)
         uL = ucL(2)/ucL(1)
         vL = ucL(3)/ucL(1)
         wL = ucL(4)/ucL(1)
        qnL = uL*nx + vL*ny + wL*nz
         pL = (gamma-one)*( ucL(5) - half*rhoL*(uL*uL+vL*vL+wL*wL) )
         aL = ddt_sqrt(gamma*pL/rhoL)
         HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL+wL*wL)
         if (ddt_sqrt(uL*uL + vL*vL + wL*wL) < epsilon * aL) then
            uR2L_ddt = (epsilon * aL)**2
         elseif (ddt_sqrt(uL*uL + vL*vL + wL*wL) < aL) then
            uR2L_ddt = uL*uL + vL*vL + wL*wL
         else
            uR2L_ddt = aL**2
         endif


   !  Right state
   
       rhoR = ucR(1)
         uR = ucR(2)/ucR(1)
         vR = ucR(3)/ucR(1)
         wR = ucR(4)/ucR(1)
        qnR = uR*nx + vR*ny + wR*nz
         pR = (gamma-one)*( ucR(5) - half*rhoR*(uR*uR+vR*vR+wR*wR) )
         aR = ddt_sqrt(gamma*pR/rhoR)
         HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR+wR*wR)
         if (ddt_sqrt(uR*uR + vR*vR + wR*wR) < epsilon * aR) then
            uR2R_ddt = (epsilon * aR)**2
         elseif (ddt_sqrt(uR*uR + vR*vR + wR*wR) < aR) then
            uR2R_ddt = uR*uR + vR*vR + wR*wR
         else
            uR2R_ddt = aR**2
         endif
         
   !Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)
   
     fL(1) = rhoL*qnL
     fL(2) = rhoL*qnL * uL + pL*nx
     fL(3) = rhoL*qnL * vL + pL*ny
     fL(4) = rhoL*qnL * wL + pL*nz
     fL(5) = rhoL*qnL * HL
   
     fR(1) = rhoR*qnR
     fR(2) = rhoR*qnR * uR + pR*nx
     fR(3) = rhoR*qnR * vR + pR*ny
     fR(4) = rhoR*qnR * wR + pR*nz
     fR(5) = rhoR*qnR * HR
   
   !First compute the Roe-averaged quantities
   
   !  NOTE: See http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
   !        the Roe-averaged density.
   
       
     rho = half * (rhoR + rhoL)                         !Arithemtic-averaged density
     u = half * (uL   + uR  )                           !Arithemtic-averaged x-velocity
     v = half * (vL   + vR  )                           !Arithemtic-averaged y-velocity
     w = half * (wL   + wR  )                           !Arithemtic-averaged z-velocity
     H = half * (HL   + hR  )                           !Arithemtic-averaged total enthalpy
     p = half * (pL   + pR  )                           !Arithmetic-averaged pressure
     a = ddt_sqrt( (gamma-one)*(H-half*(u*u + v*v + w*w)) ) !Arithemtic-averaged speed of sound
    qn = u*nx + v*ny + w*nz                             !Arithemtic-averaged face-normal velocity
  !  uR2 = half * (uR2L + uR2R)                           !Arithmetic-averaged scaling term
   uR2 = half * (uR2L_ddt + uR2R_ddt)                   !Arithmetic-averaged scaling term

   !Wave Strengths
   
      drho = rhoR - rhoL !Density difference
        dp =   pR - pL   !Pressure difference
       dqn =  qnR - qnL  !Normal velocity difference

       drhou = ucR(2) - ucL(2)
       drhov = ucR(3) - ucL(3)
       drhow = ucR(4) - ucL(4)
       drhoE = ucR(5) - ucL(5)
       absU  = ddt_abs(qn)

      ! ! Entropy fix? I guess not....... Need to investigate more
      ! if (absU < 0.1_p2 * a) then
      !     ! if ( ws(i) < dws(i) ) ws(i) = half * ( ws(i)*ws(i)/dws(i)+dws(i) )
      !     absU  = half * ( (absU**2)/(0.2_p2*a) + (0.2_p2*a) ) 
      ! end if
      
        beta = rho/(gamma*p)
       alpha = half * (one-beta*uR2)
      cprime = ddt_sqrt((alpha**2) * (qn**2) + uR2)
      uprime = qn * (one - alpha)

       cstar = half * (ddt_abs(uprime + cprime) + ddt_abs(uprime - cprime))
       Mstar = half * (ddt_abs(uprime + cprime) - ddt_abs(uprime - cprime)) / cprime

        delu = Mstar*dqn + (cstar - (one-two*alpha)*absU - alpha*qn*Mstar)*(dp/(rho*uR2))
        delp = Mstar*dp  + (cstar - absU + alpha*qn*Mstar) * rho * dqn

        diss_cartesian(:,1) = (absU * drho + delu * rho) * njk
        diss_cartesian(:,2) = (absU * drhou + delu * rho * u) * njk
        diss_cartesian(:,3) = (absU * drhov + delu * rho * v) * njk
        diss_cartesian(:,4) = (absU * drhow + delu * rho * w) * njk
        diss_cartesian(:,5) = (absU * drhoE + delu * rho * H) * njk

        ! Equivalent of matrix multiplication for the ddt type
        diss(1) = diss_cartesian(1,1)*njk(1) + diss_cartesian(2,1)*njk(2) + diss_cartesian(3,1)*njk(3)
        diss(2) = diss_cartesian(1,2)*njk(1) + diss_cartesian(2,2)*njk(2) + diss_cartesian(3,2)*njk(3)
        diss(3) = diss_cartesian(1,3)*njk(1) + diss_cartesian(2,3)*njk(2) + diss_cartesian(3,3)*njk(3)
        diss(4) = diss_cartesian(1,4)*njk(1) + diss_cartesian(2,4)*njk(2) + diss_cartesian(3,4)*njk(3)
        diss(5) = diss_cartesian(1,5)*njk(1) + diss_cartesian(2,5)*njk(2) + diss_cartesian(3,5)*njk(3)

        ! now add the final term
        ! diss(1) += 0
        diss(2) = diss(2) + delp * nx
        diss(3) = diss(3) + delp * ny
        diss(4) = diss(4) + delp * nz
        diss(5) = diss(5) + delp * qn

   ! This is the numerical flux: Roe flux = 1/2 *[  Fn(UL)+Fn(UR) - |An|(UR-UL) ]
   
     numerical_flux = half * (fL + fR - diss)
   
   !--------------
   ! Output
   
      ! Normal max wave speed
       temp = ddt_abs(qn) + a
        wsn = temp%f
   
     do i = 1, 5
   
      !Numerical flux
       num_flux(i) = numerical_flux(i)%f
   
      do j = 1, 5
   
        !Flux derivative
         dFdU(i,j) = numerical_flux(i)%df(j)
   
      end do
   
     end do
 

    !  write(*,*) 'CONS'
    !  write(*,*) fL(:)%f
    !  write(*,*) fR(:)%f
    !  write(*,*) diss(:)%f
    !  write(*,*)
    end subroutine roe_low_mach_ddt
   !--------------------------------------------------------------------------------


    subroutine roe_prim_low_mach_ddt(qcL, qcR, uR2L, uR2R, njk, num_flux,dFdU,wsn)

      use derivative_data_df5
     
      implicit none
      integer , parameter :: p2 = selected_real_kind(15) ! Double precision
     
     !Input
      type(derivative_data_type_df5), dimension(5), intent( in) :: qcL
      type(derivative_data_type_df5), dimension(5), intent( in) :: qcR
      real(p2)                                    , intent( in) :: uR2L, uR2R
      real(p2)                      , dimension(3), intent( in) :: njk
     
     !Output
      real(p2), dimension(5)  , intent(out) :: num_flux !Numerical viscous flux
      real(p2), dimension(5,5), intent(out) :: dFdU     !Numerical viscous flux Jacobian
      real(p2),                 intent(out) :: wsn      ! Max wave speed
     
     !Some constants
      real(p2) ::    zero = 0.00_p2
      real(p2) ::     one = 1.00_p2
      real(p2) ::     two = 2.00_p2
      real(p2) ::    half = 0.50_p2
      real(p2) ::   gamma = 1.40_p2                 ! Ratio of specific heats
     
     !Local variables
     !
     !            L = Left
     !            R = Right
     ! No subscript = Roe average
     
      type(derivative_data_type_df5), dimension(5) :: numerical_flux !Numerical flux in ddt
      real(p2)                       :: nx, ny, nz             ! Normal vector components
     
      type(derivative_data_type_df5) :: uL, uR, vL, vR, wL, wR ! Velocity components.
      type(derivative_data_type_df5) :: rhoL, rhoR, pL, pR     ! Primitive variables.
      type(derivative_data_type_df5) :: qnL, qnR               ! Normal velocities
      type(derivative_data_type_df5) :: aL, aR, HL, HR, EL, ER ! Speed of sound, Total enthalpy
      type(derivative_data_type_df5), dimension(5)   :: fL     ! Physical flux evaluated at ucL
      type(derivative_data_type_df5), dimension(5)   :: fR     ! Physical flux evaluated at ucR
     
      type(derivative_data_type_df5) :: RT                     ! RT = sqrt(rhoR/rhoL)
      type(derivative_data_type_df5) :: rho,u,v,w,H,a,qn,p     ! Roe-averages
     
      type(derivative_data_type_df5) :: drho, dqn, dp          ! Differences in rho, qn, p, e.g., dp=pR-pL
      type(derivative_data_type_df5) :: drhou, drhov, drhow, drhoE ! Differences in conserved vars
      type(derivative_data_type_df5) :: absU
      type(derivative_data_type_df5) :: uprime, cprime, alpha, beta
      real(p2)                       :: uR2
      type(derivative_data_type_df5) :: cstar, Mstar, delu, delp
      type(derivative_data_type_df5), dimension(4) :: LdU      ! Wave strengths = L*(UR-UL)
      type(derivative_data_type_df5) :: du, dv, dw             ! Velocity differences
      type(derivative_data_type_df5), dimension(4)   :: ws     ! Wave speeds
      type(derivative_data_type_df5), dimension(4)   :: dws    ! Width of a parabolic fit for entropy fix
      type(derivative_data_type_df5), dimension(5,4) :: R      ! Right-eigenvector matrix
      type(derivative_data_type_df5), dimension(5)   :: diss   ! Dissipation term
     
      type(derivative_data_type_df5) :: temp
      
      type(derivative_data_type_df5), dimension(3,5) :: diss_cartesian
      
      integer                        :: i, j
     
     ! Face normal vector (unit vector)
     
       nx = njk(1)
       ny = njk(2)
       nz = njk(3)
     
     !Primitive and other variables.
     
     !  Left state
     
         rhoL = qcL(1)*gamma/qcL(5)
           uL = qcL(2)
           vL = qcL(3)
           wL = qcL(4)
          qnL = uL*nx + vL*ny + wL*nz
           pL = qcL(1)
           aL = ddt_sqrt(gamma*pL/rhoL)
           HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL+wL*wL)
           EL = qcL(1)/(gamma-one)+half*rhoL*(qcL(2)*qcL(2)+qcL(3)*qcL(3)+qcL(4)*qcL(4))
     
     !  Right state
     
         rhoR = qcR(1)*gamma/qcR(5)
           uR = qcR(2)
           vR = qcR(3)
           wR = qcR(4)
          qnR = uR*nx + vR*ny + wR*nz
           pR = qcR(1)
           aR = ddt_sqrt(gamma*pR/rhoR)
           HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR+wR*wR)
           ER = qcR(1)/(gamma-one)+half*rhoR*(qcR(2)*qcR(2)+qcR(3)*qcR(3)+qcR(4)*qcR(4))
     
     !Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)
     
       fL(1) = rhoL*qnL
       fL(2) = rhoL*qnL * uL + pL*nx
       fL(3) = rhoL*qnL * vL + pL*ny
       fL(4) = rhoL*qnL * wL + pL*nz
       fL(5) = rhoL*qnL * HL
     
       fR(1) = rhoR*qnR
       fR(2) = rhoR*qnR * uR + pR*nx
       fR(3) = rhoR*qnR * vR + pR*ny
       fR(4) = rhoR*qnR * wR + pR*nz
       fR(5) = rhoR*qnR * HR
     
     !First compute the Roe-averaged quantities
     
     !  NOTE: See http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
     !        the Roe-averaged density.
     
         
       rho = half * (rhoR + rhoL)                           !Arithemtic-averaged density
       u = half * (uL   + uR  )                           !Arithemtic-averaged x-velocity
       v = half * (vL   + vR  )                           !Arithemtic-averaged y-velocity
       w = half * (wL   + wR  )                           !Arithemtic-averaged z-velocity
       H = half * (HL   + hR  )                           !Arithemtic-averaged total enthalpy
       p = half * (pL   + pR  )                           !Arithmetic-averaged pressure
       a = ddt_sqrt( (gamma-one)*(H-half*(u*u + v*v + w*w)) ) !Arithemtic-averaged speed of sound
      qn = u*nx + v*ny + w*nz                             !Arithemtic-averaged face-normal velocity
     uR2 = half * (uR2L + uR2R)                           !Arithmetic-averaged scaling term
  
     !Wave Strengths
     
        drho = rhoR - rhoL !Density difference
          dp =   pR - pL   !Pressure difference
         dqn =  qnR - qnL  !Normal velocity difference
  
         drhou = qcR(2)*rhoR - qcL(2)*rhoL
         drhov = qcR(3)*rhoR - qcL(3)*rhoL
         drhow = qcR(4)*rhoR - qcL(4)*rhoL
         drhoE = ER*rhoR - EL*rhoL
         absU  = ddt_abs(qn)
  
          beta = rho/(gamma*p)
         alpha = half * (one-beta*uR2)
        cprime = ddt_sqrt((alpha**2) * (qn**2) + uR2)
        uprime = qn * (one - alpha)
  
         cstar = half * (ddt_abs(uprime + cprime) + ddt_abs(uprime - cprime))
         Mstar = half * (ddt_abs(uprime + cprime) - ddt_abs(uprime - cprime)) / cprime
  
          delu = Mstar*dqn + (cstar - (one-two*alpha)*absU - alpha*qn*Mstar)*(dp/(rho*uR2))
          delp = Mstar*dp  + (cstar - absU + alpha*qn*Mstar) * rho * dqn
  
          diss_cartesian(:,1) = (absU * drho + delu * rho) * njk
          diss_cartesian(:,2) = (absU * drhou + delu * rho * u) * njk
          diss_cartesian(:,3) = (absU * drhov + delu * rho * v) * njk
          diss_cartesian(:,4) = (absU * drhow + delu * rho * w) * njk
          diss_cartesian(:,5) = (absU * drhoE + delu * rho * H) * njk
  
          ! Equivalent of matrix multiplication for the ddt type
          diss(1) = diss_cartesian(1,1)*njk(1) + diss_cartesian(2,1)*njk(2) + diss_cartesian(3,1)*njk(3)
          diss(2) = diss_cartesian(1,2)*njk(1) + diss_cartesian(2,2)*njk(2) + diss_cartesian(3,2)*njk(3)
          diss(3) = diss_cartesian(1,3)*njk(1) + diss_cartesian(2,3)*njk(2) + diss_cartesian(3,3)*njk(3)
          diss(4) = diss_cartesian(1,4)*njk(1) + diss_cartesian(2,4)*njk(2) + diss_cartesian(3,4)*njk(3)
          diss(5) = diss_cartesian(1,5)*njk(1) + diss_cartesian(2,5)*njk(2) + diss_cartesian(3,5)*njk(3)
  
          ! now add the final term
          ! diss(1) += 0
          diss(2) = diss(2) + delp * nx
          diss(3) = diss(3) + delp * ny
          diss(4) = diss(4) + delp * nz
          diss(5) = diss(5) + delp * qn
  
     ! This is the numerical flux: Roe flux = 1/2 *[  Fn(UL)+Fn(UR) - |An|(UR-UL) ]
     
      !  numerical_flux = half * (fL + fR - diss)
          numerical_flux = (fL + fR)
     !--------------
     ! Output
     
        ! Normal max wave speed
         temp = ddt_abs(qn) + a
          wsn = temp%f
     
       do i = 1, 5
     
        !Numerical flux
         num_flux(i) = numerical_flux(i)%f
     
        do j = 1, 5
     
          !Flux derivative
           dFdU(i,j) = numerical_flux(i)%df(j)
     
        end do
     
       end do

      !  write(*,*) 'Prim'
      !  write(*,*) fL(:)%f
      !  write(*,*) fR(:)%f
      !  write(*,*) diss(:)%f
      !  write(*,*)
      end subroutine roe_prim_low_mach_ddt
     !--------------------------------------------------------------------------------
  
!********************************************************************************
!* -- 3D HLL Flux Function and Jacobian --
!*
!*
!* This subroutine computes the HLL flux for the Euler equations
!* in the direction, njk=[nx,ny,nz].
!*
!* A. Harten, P. D. Lax, and B. van Leer,On Upstream Differencing and 
!* Godunov-Type Schemes for Hyperbolic Conservation Laws, SIAM Review,
!* 25(1), pp. 35-61, 1983.
!*
!* With wave speeds evaluated by Einfeldt's method:
!* B. Einfeldt, On Godunov-Type Methods for Gas Dynamics, SIAM Journal of
!* Numerical Analysis, 25(2), pp. 294-318, 1988.
!*
!*
!* Conservative form of the Euler equations:
!*
!*     dU/dt + dF/dx + dG/dy + dH/dz = 0
!*
!* This subroutine computes the numerical flux for the flux in the direction,
!* njk=[nx,ny,nz]:
!*
!*     Fn = F*nx + G*ny + H*nz = | rho*qn          |
!*                               | rho*qn*u + p*nx |
!*                               | rho*qn*v + p*ny |
!*                               | rho*qn*w + p*nz |
!*                               | rho*qn*H        |    (qn = u*nx + v*ny + w*nz)
!*
!* The HLL flux is implemented in the following form:
!*
!*   Numerical flux = SRp*Fn(UR)-SLm*Fn(UL) + SLm*SRp*(UR-UL))/(SRp-SLm)
!*
!* where the wave speeds are evaluated by Einfeldt's method,
!*
!*  SLm = min( (qn - a)_L, (qn - a)_Roe, 0)
!*  SRp = max( (qn + a)_R, (qn + a)_Roe, 0)
!*
!* Note: In the code, the vector of conserative variables are denoted by uc.
!*
!* ------------------------------------------------------------------------------
!*  Input: ucL(1:5) =  Left state (rhoL, rhoL*uL, rhoL*vL, rhoL*wR, rhoL*EL)
!*         ucR(1:5) = Right state (rhoR, rhoL*uR, rhoL*vR, rhoL*wR, rhoL*ER)
!*         njk(1:3) = unit face normal vector (nx, ny, nz), pointing from Left to Right.
!*
!*           njk
!*  Face normal ^   o Right data point
!*              |  .
!*              | .
!*              |. 
!*       -------x-------- Face
!*             .                 Left and right states are
!*            .                   1. Values at data points for 1st-order accuracy
!*           .                    2. Extrapolated values at the face midpoint 'x'
!*          o Left data point        for 2nd/higher-order accuracy.
!*
!*
!* Output:  num_flux(1:5) = the numerical flux vector
!*          dFdU(1:5,1:5) = flux Jacobian matrix
!*                    wsn = maximum wave speed (eigenvalue)
!*
!* ------------------------------------------------------------------------------
!*
!* Note: This function is a ddt-version, which means that each variable carries
!*       its derivatives, and that the resulting flux "numerical_flux" will have
!*       its derivatives in "numerical_flux%df".
!*
!* Note: This subroutine has been prepared for an educational purpose.
!*       It is not at all efficient. Think about how you can optimize it.
!*       One way to make it efficient is to reduce the number of local variables,
!*       by re-using temporary variables as many times as possible.
!*
!* Note: Please let me know if you find bugs. I'll greatly appreciate it and
!*       fix the bugs.
!*
!* Katate Masatsuka, November 2012. http://www.cfdbooks.com
!********************************************************************************
 subroutine hll_ddt(ucL, ucR, njk, num_flux,dFdU,wsn)

 use derivative_data_df5

 implicit none

 integer , parameter :: p2 = selected_real_kind(15) ! Double precision

!Input
 type(derivative_data_type_df5), dimension(5), intent( in) :: ucL
 type(derivative_data_type_df5), dimension(5), intent( in) :: ucR
 real(p2)                      , dimension(3), intent( in) :: njk

!Output
 real(p2), dimension(5)  , intent(out) :: num_flux !Numerical viscous flux
 real(p2), dimension(5,5), intent(out) :: dFdU     !Numerical viscous flux Jacobian
 real(p2),                 intent(out) :: wsn      ! Max wave speed


!Some constants
 real(p2) ::   zero = 0.00_p2
 real(p2) ::    one = 1.00_p2
 real(p2) ::   half = 0.50_p2
 real(p2) ::  gamma = 1.40_p2                 ! Ratio of specific heats

!Local variables
 real(p2)                       :: nx, ny, nz             ! Normal vector components
 type(derivative_data_type_df5), dimension(5) :: numerical_flux
 type(derivative_data_type_df5) :: uL, uR, vL, vR, wL, wR ! Velocity components.
 type(derivative_data_type_df5) :: rhoL, rhoR, pL, pR     ! Primitive variables.
 type(derivative_data_type_df5) :: qnL, qnR               ! Normal velocities
 type(derivative_data_type_df5) :: aL, aR, HL, HR         ! Speed of sound, Total enthalpy
 type(derivative_data_type_df5) :: RT,rho,u,v,w,H,a,qn    ! Roe-averages
 type(derivative_data_type_df5) :: SLm, SRp               ! Einfeldt's wave speed estimates
 type(derivative_data_type_df5), dimension(5)   :: fL     ! Physical flux evaluated at ucL
 type(derivative_data_type_df5), dimension(5)   :: fR     ! Physical flux evaluated at ucR

 type(derivative_data_type_df5) :: temp
 integer                        :: i, j

! Face normal vector (unit vector)

  nx = njk(1)
  ny = njk(2)
  nz = njk(3)

!Primitive and other variables.

!  Left state

    rhoL = ucL(1)
      uL = ucL(2)/ucL(1)
      vL = ucL(3)/ucL(1)
      wL = ucL(4)/ucL(1)
     qnL = uL*nx + vL*ny + wL*nz
      pL = (gamma-one)*( ucL(5) - half*rhoL*(uL*uL+vL*vL+wL*wL) )
      aL = ddt_sqrt(gamma*pL/rhoL)
      HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL+wL*wL)

!  Right state

    rhoR = ucR(1)
      uR = ucR(2)/ucR(1)
      vR = ucR(3)/ucR(1)
      wR = ucR(4)/ucR(1)
     qnR = uR*nx + vR*ny + wR*nz
      pR = (gamma-one)*( ucR(5) - half*rhoR*(uR*uR+vR*vR+wR*wR) )
      aR = ddt_sqrt(gamma*pR/rhoR)
      HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR+wR*wR)

! Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)

  fL(1) = rhoL*qnL
  fL(2) = rhoL*qnL * uL + pL*nx
  fL(3) = rhoL*qnL * vL + pL*ny
  fL(4) = rhoL*qnL * wL + pL*nz
  fL(5) = rhoL*qnL * HL

  fR(1) = rhoR*qnR
  fR(2) = rhoR*qnR * uR + pR*nx
  fR(3) = rhoR*qnR * vR + pR*ny
  fR(4) = rhoR*qnR * wR + pR*nz
  fR(5) = rhoR*qnR * HR

! Compute the Roe-averaged quantities
!  NOTE: See http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
!        the Roe-averaged density.

    RT = ddt_sqrt(rhoR/rhoL)
   rho = RT*rhoL                                            !Roe-averaged density
     u = (uL + RT*uR)/(one + RT)                            !Roe-averaged x-velocity
     v = (vL + RT*vR)/(one + RT)                            !Roe-averaged y-velocity
     w = (wL + RT*wR)/(one + RT)                            !Roe-averaged z-velocity
     H = (HL + RT*HR)/(one + RT)                            !Roe-averaged total enthalpy
     a = ddt_sqrt( (gamma-one)*(H-half*(u*u + v*v + w*w)) ) !Roe-averaged speed of sound
    qn = u*nx + v*ny + w*nz                                 !Roe-averaged face-normal velocity

! Einfeldt's wave speed estimates

   SLm = ddt_min( ddt_min(qnL-aL, qn-a), zero)
   SRp = ddt_max( ddt_max(qnR+aR, qn+a), zero)

! Compute the HLL flux

  numerical_flux = (SRp*fL-SLm*fR + SLm*SRp*(ucR-ucL))/(SRp-SLm)

!--------------
! Output

   ! Normal max wave speed
    temp = ddt_abs(qn) + a
     wsn = temp%f

  do i = 1, 5

   !Numerical flux
    num_flux(i) = numerical_flux(i)%f

   do j = 1, 5

     !Flux derivative
      dFdU(i,j) = numerical_flux(i)%df(j)

   end do

  end do

 end subroutine hll_ddt
!--------------------------------------------------------------------------------

!********************************************************************************
!* -- 3D Rusanov Flux Function and Jacobian --
!*
!*
!* This subroutine computes the Roe flux for the Euler equations
!* in the direction, njk=[nx,ny,nz].
!*
!* V. V. Rusanov, Calculation of Interaction of Non-Steady Shock Waves with
!* Obstacles, J. Comput. Math. Phys. USSR, 1, pp. 267-279, 1961.
!*
!*
!* Conservative form of the Euler equations:
!*
!*     dU/dt + dF/dx + dG/dy + dH/dz = 0
!*
!* This subroutine computes the numerical flux for the flux in the direction,
!* njk=[nx,ny,nz]:
!*
!*     Fn = F*nx + G*ny + H*nz = | rho*qn          |
!*                               | rho*qn*u + p*nx |
!*                               | rho*qn*v + p*ny |
!*                               | rho*qn*w + p*nz |
!*                               | rho*qn*H        |    (qn = u*nx + v*ny + w*nz)
!*
!* The Rusanov flux is implemented in the following form:
!*
!*   Numerical flux = 1/2 [ Fn(UR) + Fn(UL) - max_eigenvalue*dU ]
!*
!* where
!*
!*  max_eigenvalue = |qn| + a,
!*
!* evaluated at the Roe-average state.
!*
!* Note: In the code, the vector of conserative variables are denoted by uc.
!*
!* ------------------------------------------------------------------------------
!*  Input: ucL(1:5) =  Left state (rhoL, rhoL*uL, rhoL*vL, rhoL*wR, rhoL*EL)
!*         ucR(1:5) = Right state (rhoR, rhoL*uR, rhoL*vR, rhoL*wR, rhoL*ER)
!*         njk(1:3) = unit face normal vector (nx, ny, nz), pointing from Left to Right.
!*
!*           njk
!*  Face normal ^   o Right data point
!*              |  .
!*              | .
!*              |. 
!*       -------x-------- Face
!*             .                 Left and right states are
!*            .                   1. Values at data points for 1st-order accuracy
!*           .                    2. Extrapolated values at the face midpoint 'x'
!*          o Left data point        for 2nd/higher-order accuracy.
!*
!*
!* Output:  num_flux(1:5) = the numerical flux vector
!*          dFdU(1:5,1:5) = flux Jacobian matrix
!*                    wsn = maximum wave speed (eigenvalue)
!*
!* ------------------------------------------------------------------------------
!*
!* Note: This function is a ddt-version, which means that each variable carries
!*       its derivatives, and that the resulting flux "numerical_flux" will have
!*       its derivatives in "numerical_flux%df".
!*
!* Note: This subroutine has been prepared for an educational purpose.
!*       It is not at all efficient. Think about how you can optimize it.
!*       One way to make it efficient is to reduce the number of local variables,
!*       by re-using temporary variables as many times as possible.
!*
!* Note: Please let me know if you find bugs. I'll greatly appreciate it and
!*       fix the bugs.
!*
!* Katate Masatsuka, November 2012. http://www.cfdbooks.com
!********************************************************************************
 subroutine rusanov_ddt(ucL, ucR, njk, num_flux,dFdU,wsn)

 use derivative_data_df5

 implicit none

 integer , parameter :: p2 = selected_real_kind(15) ! Double precision

!Input
 type(derivative_data_type_df5), dimension(5), intent( in) :: ucL
 type(derivative_data_type_df5), dimension(5), intent( in) :: ucR
 real(p2)                      , dimension(3), intent( in) :: njk

!Output
 real(p2), dimension(5)  , intent(out) :: num_flux !Numerical viscous flux
 real(p2), dimension(5,5), intent(out) :: dFdU     !Numerical viscous flux Jacobian
 real(p2),                 intent(out) :: wsn      ! Max wave speed

!Some constants
 real(p2) ::    one = 1.00_p2
 real(p2) ::   half = 0.50_p2
 real(p2) ::  gamma = 1.40_p2                 ! Ratio of specific heats

!Local variables
 real(p2)                       :: nx, ny, nz             ! Normal vector components
 type(derivative_data_type_df5), dimension(5) :: numerical_flux !Numerical viscous flux
 type(derivative_data_type_df5) :: uL, uR, vL, vR, wL, wR ! Velocity components.
 type(derivative_data_type_df5) :: rhoL, rhoR, pL, pR     ! Primitive variables.
 type(derivative_data_type_df5) :: qnL, qnR               ! Normal velocities
 type(derivative_data_type_df5) :: aL, aR, HL, HR         ! Speed of sound, Total enthalpy
 type(derivative_data_type_df5) :: RT,rho,u,v,w,H,a,qn    ! Roe-averages
 type(derivative_data_type_df5) :: max_eigenvalue         ! Maximum eigenvalue = (|qn| + a) at Roe averaged state
 type(derivative_data_type_df5), dimension(5)   :: fL     ! Physical flux evaluated at ucL
 type(derivative_data_type_df5), dimension(5)   :: fR     ! Physical flux evaluated at ucR

 type(derivative_data_type_df5) :: temp
 integer                        :: i, j

! Face normal vector (unit vector)

  nx = njk(1)
  ny = njk(2)
  nz = njk(3)

!Primitive and other variables.

!  Left state

    rhoL = ucL(1)
      uL = ucL(2)/ucL(1)
      vL = ucL(3)/ucL(1)
      wL = ucL(4)/ucL(1)
     qnL = uL*nx + vL*ny + wL*nz
      pL = (gamma-one)*( ucL(5) - half*rhoL*(uL*uL+vL*vL+wL*wL) )
      aL = ddt_sqrt(gamma*pL/rhoL)
      HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL+wL*wL)

!  Right state

    rhoR = ucR(1)
      uR = ucR(2)/ucR(1)
      vR = ucR(3)/ucR(1)
      wR = ucR(4)/ucR(1)
     qnR = uR*nx + vR*ny + wR*nz
      pR = (gamma-one)*( ucR(5) - half*rhoR*(uR*uR+vR*vR+wR*wR) )
      aR = ddt_sqrt(gamma*pR/rhoR)
      HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR+wR*wR)

! Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)

  fL(1) = rhoL*qnL
  fL(2) = rhoL*qnL * uL + pL*nx
  fL(3) = rhoL*qnL * vL + pL*ny
  fL(4) = rhoL*qnL * wL + pL*nz
  fL(5) = rhoL*qnL * HL

  fR(1) = rhoR*qnR
  fR(2) = rhoR*qnR * uR + pR*nx
  fR(3) = rhoR*qnR * vR + pR*ny
  fR(4) = rhoR*qnR * wR + pR*nz
  fR(5) = rhoR*qnR * HR

! Compute the Roe-averaged quantities
!  NOTE: See http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
!        the Roe-averaged density.

    RT = ddt_sqrt(rhoR/rhoL)
   rho = RT*rhoL                                            !Roe-averaged density
     u = (uL + RT*uR)/(one + RT)                            !Roe-averaged x-velocity
     v = (vL + RT*vR)/(one + RT)                            !Roe-averaged y-velocity
     w = (wL + RT*wR)/(one + RT)                            !Roe-averaged z-velocity
     H = (HL + RT*HR)/(one + RT)                            !Roe-averaged total enthalpy
     a = ddt_sqrt( (gamma-one)*(H-half*(u*u + v*v + w*w)) ) !Roe-averaged speed of sound
    qn = u*nx + v*ny + w*nz                                 !Roe-averaged face-normal velocity

! Maximum of absolute values of eigenvalues:

  max_eigenvalue = ( ddt_abs(qn) + a )

! Local Lax-Friedrichs flux = 1/2 *[  Fn(UL)+Fn(UR) - |maximum eigenvalue|(UR-UL) ]

  numerical_flux = half * (fL + fR - max_eigenvalue*(ucR- ucL) )

!--------------
! Output

   ! Normal max wave speed
    temp = ddt_abs(qn) + a
     wsn = temp%f

  do i = 1, 5

   !Numerical flux
    num_flux(i) = numerical_flux(i)%f

   do j = 1, 5

     !Flux derivative
      dFdU(i,j) = numerical_flux(i)%df(j)

   end do

  end do

 end subroutine rusanov_ddt
!--------------------------------------------------------------------------------

!********************************************************************************
!* -- 3D Rotated-Roe-HLL (Rotated-RHLL) Flux Function and Jacobian ---
!*
!* This subroutine computes the Rotated-RHLL flux for the Euler equations
!* in the direction, njk=[nx,ny,nz].
!*
!*
!* H. Nishikawa and K. Kitamura, Very Simple, Carbuncle-Free, Boundary-Layer
!* Resolving, Rotated-Hybrid Riemann Solvers,
!* Journal of Computational Physics, 227, pp. 2560-2581, 2008.
!* https://doi.org/10.1016/j.jcp.2007.11.003
!* http://ossanworld.com/hiroakinishikawa/My_papers/nishikawa_kitamura_jcp2008v227pp2560-2581_preprint.pdf
!*
!*
!*
!* Robust for nonlinear instability (carbuncle).
!* Recommended for high-speed flows involving strong shocks.
!* It is also capable of resolving boundary layers.
!*
!*
!* Conservative form of the Euler equations:
!*
!*     dU/dt + dF/dx + dG/dy + dH/dz = 0
!*
!* This subroutine computes the numerical flux for the flux in the direction,
!* njk=[nx,ny,nz]:
!*
!*     Fn = F*nx + G*ny + H*nz = | rho*qn          |
!*                               | rho*qn*u + p*nx |
!*                               | rho*qn*v + p*ny |
!*                               | rho*qn*w + p*nz |
!*                               | rho*qn*H        |    (qn = u*nx + v*ny + w*nz)
!*
!* The Rotated-RHLL flux is defined as
!*
!*   Numerical flux = alpha1*HLL(n1) + alpha2*Roe(n2),  if |dq| > eps
!*                    Roe(n)                         ,  if |dq| < eps
!*
!* where n1 is taken as the velocity difference vector (normal to shock or
!* tangent to shear waves), n2 is a vector perpendicular to n1,
!* alpha1 = n*n1, alpha2 = n*n2. That is, we decompose the normal vector, n,
!* in the two directions, n1 and n2, and apply the HLL in n1 and the Roe in n2.
!* However, the resulting flux can be simplified and can be made to look like
!* a single Roe-type flux. The additional cost over the Roe flux is reported
!* to be as small as 14% of the Roe flux.
!*
!* Note: HLL is introduced only near discontinuities and only by a fraction, alpha1.
!*
!* Note: For full convergence to machine zero for staedy state computations,
!*       the vectors, n1 and n2, may have to be frozen (by storing alpha1 and alpha2).
!*       For time-accurate calculation, it would not be necessary.
!*
!* Note: Here, the Roe flux is computed without tangent vectors.
!*
!* ------------------------------------------------------------------------------
!*  Input:   ucL(1:5) =  Left state (rhoL, rhoL*uL, rhoL*vL, rhoL*wL, rhoL*EL)
!*           ucR(1:5) = Right state (rhoR, rhoR*uR, rhoR*vR, rhoR*wR, rhoR*ER)
!*           njk(1:3) = unit face normal vector (nx, ny, nz), pointing from Left to Right
!*
!*           njk
!*  Face normal ^   o Right data point
!*              |  .
!*              | .
!*              |. 
!*       -------x-------- Face
!*             .                 Left and right states are
!*            .                   1. Values at data points for 1st-order accuracy
!*           .                    2. Extrapolated values at the face midpoint 'x'
!*          o Left data point        for 2nd/higher-order accuracy.
!*
!* Output:  num_flux(1:5) = the numerical flux vector
!*          dFdU(1:5,1:5) = flux Jacobian matrix
!*                    wsn = maximum wave speed (eigenvalue)
!*              exact_jac = .false. for non-exact Jac, = .true. for exact Jac
!*
!* ------------------------------------------------------------------------------
!*
!* Note: This function is a ddt-version, which means that each variable carries
!*       its derivatives, and that the resulting flux "numerical_flux" will have
!*       its derivatives in "numerical_flux%df". Jacobian is EXACT.
!*
!* Note: This subroutine has been prepared for an educational purpose.
!*       It is not at all efficient. Think about how you can optimize it.
!*       One way to make it efficient is to reduce the number of local variables,
!*       by re-using temporary variables as many times as possible.
!*
!* Note: Please let me know if you find bugs. I'll greatly appreciate it and
!*       fix the bugs.
!* 
!* Katate Masatsuka, December 2012. http://www.cfdbooks.com
!********************************************************************************
 subroutine rhll_ddt(ucL, ucR, njk, num_flux,dFdU,wsn,exact_jac)

 use derivative_data_df5

 implicit none
 integer , parameter :: p2 = selected_real_kind(15) ! Double precision

!Input
 type(derivative_data_type_df5), dimension(5), intent( in) :: ucL
 type(derivative_data_type_df5), dimension(5), intent( in) :: ucR
 real(p2)                      , dimension(3), intent( in) :: njk
 logical                                     , intent( in) :: exact_jac

!Output
 real(p2), dimension(5)  , intent(out) :: num_flux !Numerical viscous flux
 real(p2), dimension(5,5), intent(out) :: dFdU     !Numerical viscous flux Jacobian
 real(p2),                 intent(out) :: wsn      ! Max wave speed

!Some constants
 real(p2) ::    zero = 0.00_p2
 real(p2) ::     one = 1.00_p2
 real(p2) ::     two = 2.00_p2
 real(p2) ::    half = 0.50_p2
 real(p2) ::   gamma = 1.40_p2                 ! Ratio of specific heats

!Local variables
!
!            L = Left
!            R = Right
! No subscript = Roe average


 type(derivative_data_type_df5), dimension(5) :: numerical_flux !Numerical flux in ddt
 real(p2)                       :: nx, ny, nz             ! Normal vector components

 type(derivative_data_type_df5) :: uL, uR, vL, vR, wL, wR ! Velocity components.
 type(derivative_data_type_df5) :: rhoL, rhoR, pL, pR     ! Primitive variables.
 type(derivative_data_type_df5) :: qnL, qnR               ! Normal velocities
 type(derivative_data_type_df5) :: aL, aR, HL, HR         ! Speed of sound, Total enthalpy
 type(derivative_data_type_df5), dimension(5)   :: fL     ! Physical flux evaluated at ucL
 type(derivative_data_type_df5), dimension(5)   :: fR     ! Physical flux evaluated at ucR

 type(derivative_data_type_df5) :: RT                     ! RT = sqrt(rhoR/rhoL)
 type(derivative_data_type_df5) :: rho,u,v,w,H,a,qn       ! Roe-averages

 type(derivative_data_type_df5) :: drho, dqn, dp          ! Differences in rho, qn, p, e.g., dp=pR-pL
 type(derivative_data_type_df5), dimension(4) :: LdU      ! Wave strengths = L*(UR-UL)
 type(derivative_data_type_df5) :: du, dv, dw             ! Velocity differences
 type(derivative_data_type_df5), dimension(4)   :: ws     ! Wave speeds
 type(derivative_data_type_df5), dimension(4)   :: dws    ! Width of a parabolic fit for entropy fix
 type(derivative_data_type_df5), dimension(5,4) :: R      ! Right-eigenvector matrix
 type(derivative_data_type_df5), dimension(5)   :: diss   ! Dissipation term

 type(derivative_data_type_df5), dimension(4)   :: eig    ! Wave speeds (eigenvalues)

 type(derivative_data_type_df5) :: SRp,SLm                   ! Wave speeds for the HLL part
 type(derivative_data_type_df5) :: nx1, ny1, nz1             ! Vector along which HLL is applied
 type(derivative_data_type_df5) :: nx2, ny2, nz2             ! Vector along which Roe is applied
 type(derivative_data_type_df5) :: alpha1, alpha2            ! Projections of the new normals
 type(derivative_data_type_df5) :: abs_dq                    ! Magnitude of the velocity difference
 type(derivative_data_type_df5) :: temp                      ! Temporary variable

 integer                        :: i, j

! Face normal vector (unit vector)

  nx = njk(1)
  ny = njk(2)
  nz = njk(3)

!Primitive and other variables.

!  Left state

    rhoL = ucL(1)
      uL = ucL(2)/ucL(1)
      vL = ucL(3)/ucL(1)
      wL = ucL(4)/ucL(1)
     qnL = uL*nx + vL*ny + wL*nz
      pL = (gamma-one)*( ucL(5) - half*rhoL*(uL*uL+vL*vL+wL*wL) )
      aL = ddt_sqrt(gamma*pL/rhoL)
      HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL+wL*wL)

!  Right state

    rhoR = ucR(1)
      uR = ucR(2)/ucR(1)
      vR = ucR(3)/ucR(1)
      wR = ucR(4)/ucR(1)
     qnR = uR*nx + vR*ny + wR*nz
      pR = (gamma-one)*( ucR(5) - half*rhoR*(uR*uR+vR*vR+wR*wR) )
      aR = ddt_sqrt(gamma*pR/rhoR)
      HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR+wR*wR)

!Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)

  fL(1) = rhoL*qnL
  fL(2) = rhoL*qnL * uL + pL*nx
  fL(3) = rhoL*qnL * vL + pL*ny
  fL(4) = rhoL*qnL * wL + pL*nz
  fL(5) = rhoL*qnL * HL

  fR(1) = rhoR*qnR
  fR(2) = rhoR*qnR * uR + pR*nx
  fR(3) = rhoR*qnR * vR + pR*ny
  fR(4) = rhoR*qnR * wR + pR*nz
  fR(5) = rhoR*qnR * HR

!--------------------------------------------------------------------------------
! Define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
!
    abs_dq = ddt_sqrt( (uR-uL)**2 + (vR-vL)**2 + (wR-wL)**2 )

  n1_dq : if ( abs_dq%f > 1.0e-12_p2 ) then

  !----------------------------------------------------------------------
  ! n1 = Velocity difference vector: normal to shock or tangent to shear
  !----------------------------------------------------------------------

   ! n1(HLL) = unit velocity-difference vector

         nx1 = (uR-uL)/abs_dq
         ny1 = (vR-vL)/abs_dq
         nz1 = (wR-wL)/abs_dq
      alpha1 = nx*nx1 + ny*ny1 + nz*nz1 ! alpha1 = n*n1

     ! Make alpha1 always positive.
        if ( alpha1%f < zero ) then
            nx1 = -nx1
            ny1 = -ny1
            nz1 = -nz1
         alpha1 = -alpha1
        endif

   ! n2(Roe) = unit vector perpendicular to n1.
   !         = sum of (-nz1,0,nx), (0,-nz,ny), (-ny1,nx1,0).

         temp = ddt_sqrt( (-nz1-ny1)**2 + (-nz1+nx1)**2 + (ny1+nx1)**2 )
          nx2 = (-nz1      -ny1)/temp
          ny2 = (     -nz1 +nx1)/temp
          nz2 = ( nx1 +ny1     )/temp
       alpha2 = nx*nx2 + ny*ny2 + nz*nz2 ! alpha2 = n*n2

      ! Make alpha2 always positive.
        if ( alpha2%f < zero ) then
            nx2 = -nx2
            ny2 = -ny2
            nz2 = -nz2
         alpha2 = -alpha2
        endif

  else

  !----------------------------------------------------------------------
  ! n1 = Face tangent vector: RHLL = Roe in this case.
  !----------------------------------------------------------------------

   ! n1(HLL) = Face tangent vector (no contribution)

        temp = sqrt( (-nz - ny)**2 + (-nz + nx)**2 + (ny + nx)**2 )
         nx1 = (-nz - ny)/temp
         ny1 = (-nz + nx)/temp
         nz1 = ( ny + nx)/temp
      alpha1 = zero ! alpha1 = n*n1 = 0

   ! n2(Roe) = n
         nx2 = nx
         ny2 = ny
         nz2 = nz
      alpha2 = one  ! alpha2 = n*n2 = 1

  endif n1_dq

!--------------------------------------------------------------------------------

 ! Turn off the derivatives of these parameters if "exact_jac" = .false.
 ! This makes the Jacobian a combination of
 ! exact Roe and HLL Jacobians (approximate Jacobian for RHLL).

  if (.not.exact_jac) then
       nx1%df = zero
       ny1%df = zero
       nz1%df = zero
    alpha1%df = zero

       nx2%df = zero
       ny2%df = zero
       nz2%df = zero
    alpha2%df = zero
  endif

!--------------------------------------------------------------------------------
! Now we are going to compute the Roe flux with n2 as the normal with modified 
! wave speeds (5.12). NOTE: the Roe flux here is computed without tangent vectors.
! See "I do like CFD, VOL.1" for details: page 57, Equation (3.6.31).

! First compute the Roe-averaged quantities

!  NOTE: See http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
!        the Roe-averaged density.

    RT = ddt_sqrt(rhoR/rhoL)
   rho = RT*rhoL                                            !Roe-averaged density
     u = (uL + RT*uR)/(one + RT)                            !Roe-averaged x-velocity
     v = (vL + RT*vR)/(one + RT)                            !Roe-averaged y-velocity
     w = (wL + RT*wR)/(one + RT)                            !Roe-averaged z-velocity
     H = (HL + RT*HR)/(one + RT)                            !Roe-averaged total enthalpy
     a = ddt_sqrt( (gamma-one)*(H-half*(u*u + v*v + w*w)) ) !Roe-averaged speed of sound

!----------------------------------------------------
! Compute the wave speed estimates for the HLL part,
! following Einfeldt:
!
! B. Einfeldt, On Godunov-type methods for gas dynamics,
! SIAM Journal on Numerical Analysis 25 (2) (1988) 294-318.
!
! Note: HLL is applied to n1, which requires only the
!       following. See JCP2008 paper.
!
! Note: This flux is identical to HLL if alpha2=0, i.e.,
!       if n1 = n.
!

     qn  = u *nx1 + v *ny1 + w *nz1
     qnL = uL*nx1 + vL*ny1 + wL*nz1
     qnR = uR*nx1 + vR*ny1 + wR*nz1
     SLm = ddt_min( ddt_min(qnL-aL, qn-a), zero) !Minimum wave speed estimate
     SRp = ddt_max( ddt_max(qnR+aR, qn+a), zero) !Maximum wave speed estimate

! This is the only place where n1=(nx1,ny1,nz1) is used.
! n1=(nx1,ny1,nz1) is never used below.
!----------------------------------------------------

!----------------------------------------------------
! Compute the Roe flux dissipation based on n2, with
! the RHLL wave speeds.

! Wave Strengths

     qn  = u *nx2 + v *ny2 + w *nz2
     qnL = uL*nx2 + vL*ny2 + wL*nz2
     qnR = uR*nx2 + vR*ny2 + wR*nz2

   drho = rhoR - rhoL !Density difference
     dp =   pR - pL   !Pressure difference
    dqn =  qnR - qnL  !Normal velocity difference

  LdU(1) = (dp - rho*a*dqn )/(two*a*a) !Left-moving acoustic wave strength
  LdU(2) =  drho - dp/(a*a)            !Entropy wave strength
  LdU(3) = (dp + rho*a*dqn )/(two*a*a) !Right-moving acoustic wave strength
  LdU(4) = rho                         !Shear wave strength (not really, just a factor)

! Wave Speed (Eigenvalues)

  eig(1) = qn-a !Left-moving acoustic wave velocity
  eig(2) = qn   !Entropy wave velocity
  eig(3) = qn+a !Right-moving acoustic wave velocity
  eig(4) = qn   !Shear wave velocity

! Absolute values of the wave speeds (absolute eigenvalues)

  ws(1) = ddt_abs(qn-a) !Left-moving acoustic wave
  ws(2) = ddt_abs(qn)   !Entropy wave
  ws(3) = ddt_abs(qn+a) !Right-moving acoustic wave
  ws(4) = ddt_abs(qn)   !Shear waves

! Harten's Entropy Fix JCP(1983), 49, pp357-393. This is typically applied
! only for the nonlinear fields (k=1 and 3), but here it is applied to all
! for robustness, avoiding vanishing wave speeds by making a parabolic fit
! near ws = 0 for all waves.

  do i = 1, 4
   dws(i) = 0.2_p2*a
    if ( ws(i) < dws(i) ) ws(i) = half * ( ws(i)*ws(i)/dws(i)+dws(i) )
  end do

! Combine the wave speeds for Rotated-RHLL: Eq.(5.12) in the original JCP2008 paper.

      ws = alpha2*ws - (alpha1*two*SRp*SLm + alpha2*(SRp+SLm)*eig)/(SRp-SLm)

! Below, we compute the Roe dissipation term in the direction n2
! with the above modified wave speeds. HLL wave speeds act something like
! the entropy fix or eigenvalue limiting; they contribute only by the amount
! given by the fraction, alpha1 (less than or equal to 1.0). See JCP2008 paper.

! Right Eigenvectors:
! Note: Two shear wave components are combined into one, so that tangent vectors
!       are not required. And that's why there are only 4 vectors here.

! Left-moving acoustic wave
  R(1,1) = one    
  R(2,1) = u - a*nx2
  R(3,1) = v - a*ny2
  R(4,1) = w - a*nz2
  R(5,1) = H - a*qn

! Entropy wave
  R(1,2) = one
  R(2,2) = u
  R(3,2) = v 
  R(4,2) = w
  R(5,2) = half*(u*u + v*v + w*w)

! Right-moving acoustic wave
  R(1,3) = one
  R(2,3) = u + a*nx2
  R(3,3) = v + a*ny2
  R(4,3) = w + a*nz2
  R(5,3) = H + a*qn

! Two shear wave components combined into one (wave strength incorporated).
  du = uR - uL
  dv = vR - vL
  dw = wR - wL
  R(1,4) = zero
  R(2,4) = du - dqn*nx2
  R(3,4) = dv - dqn*ny2
  R(4,4) = dw - dqn*nz2
  R(5,4) = u*du + v*dv + w*dw - qn*dqn

!Dissipation Term: Roe dissipation with the modified wave speeds.
! |An|dU = R|Lambda|L*dU = sum_k of [ ws(k) * R(:,k) * L*dU(k) ], where n=n2.

 diss(:) = ws(1)*LdU(1)*R(:,1) + ws(2)*LdU(2)*R(:,2) &
         + ws(3)*LdU(3)*R(:,3) + ws(4)*LdU(4)*R(:,4)

!Compute the RHLL flux. (It looks like the HLL flux with Roe dissipation.)

  numerical_flux = (SRp*fL - SLm*fR)/(SRp-SLm) - half*diss

!--------------
! Output

   ! Normal max wave speed
     temp = ddt_abs(u*nx + v*ny + w*nz) + a
      wsn = temp%f

  do i = 1, 5

   !Numerical flux
    num_flux(i) = numerical_flux(i)%f

   do j = 1, 5

     !Flux derivative
      dFdU(i,j) = numerical_flux(i)%df(j)

   end do

  end do

 end subroutine rhll_ddt
!--------------------------------------------------------------------------------

 subroutine viscous_alpha_ddt(ucL,ucR,gradwL,gradwR, njk,ejk,mag_ejk, numerical_flux)

  use derivative_data_df5
  use module_input_parameter , only : M_inf, Reynolds, C_0, Freestream_Temp

  implicit none
  integer , parameter :: p2 = selected_real_kind(15)
 
 !Input
  type(derivative_data_type_df5), dimension(5)  , intent( in) :: ucL     !Left state (conservative)
  type(derivative_data_type_df5), dimension(5)  , intent( in) :: ucR     !Right state (conservative)
  real(p2)                      , dimension(3,5), intent( in) :: gradwL  !Left gradient (primitive)
  real(p2)                      , dimension(3,5), intent( in) :: gradwR  !Right gradient (primitive)
  real(p2)                      , dimension(3)  , intent( in) :: njk     !Unit directed area vector
  real(p2)                      , dimension(3)  , intent( in) :: ejk     !Unit edge vector
  real(p2)                      ,                 intent( in) :: mag_ejk !Magnitude of the edge vector
 
 !Output
  type(derivative_data_type_df5), dimension(5), intent(out) :: numerical_flux !Numerical viscous flux
 
 !Some constants
  real(p2) ::       zero = 0.0_p2
  real(p2) ::        one = 1.0_p2
  real(p2) ::      three = 3.0_p2
  real(p2) ::       half = 0.5_p2
  real(p2) ::  two_third = 2.0_p2/3.0_p2
  real(p2) :: four_third = 4.0_p2/3.0_p2
 
  real(p2) ::          C! = C_0   !Parameter for Sutherland's law
  real(p2) ::      gamma = 1.4_p2     !Ratio of specific heats
  real(p2) ::    Prandtl = 0.72_p2    !Prandtl number
  real(p2) ::     Re_inf! = Reynolds   !Free stream Reynolds number
  real(p2) ::      T_inf! = Freestream_Temp   !Free stream temperature
  ! real(p2) ::      M_inf = 2.0_p2     !Free stream Mach number
 
 !Local variables
  real(p2) :: alpha !Damping coefficient (see Nishikawa AIAA2010-5093)
  real(p2) :: Lr    !Length scale        (see Nishikawa AIAA2010-5093)
 
  type(derivative_data_type_df5) ::   uL, uR            ! x-velocity  (Left and Right states)
  type(derivative_data_type_df5) ::   vL, vR            ! y-velocity  (Left and Right states)
  type(derivative_data_type_df5) ::   wL, wR            ! z-velocity  (Left and Right states)
  type(derivative_data_type_df5) :: rhoL, rhoR          ! Density     (Left and Right states)
  type(derivative_data_type_df5) ::   pL, pR            ! Pressure    (Left and Right states)
  type(derivative_data_type_df5) ::   TL, TR            ! Temperature (Left and Right states)
 
  type(derivative_data_type_df5) :: tauxx, tauyy, tauzz !Viscous stresses: diagonal compontens
  type(derivative_data_type_df5) :: tauxy, tauyz, tauzx !Viscous stresses: off-diagonal components
  type(derivative_data_type_df5) :: tauyx, tauzy, tauxz !Viscous stresses: same as above by symmetry
  type(derivative_data_type_df5) :: qx, qy, qz          !Heat flux components
  type(derivative_data_type_df5) :: tauxn, tauyn, tauzn !Normal stresses
  type(derivative_data_type_df5) :: qn                  !Normal heat flux
 
  type(derivative_data_type_df5) :: u, v, w, T, mu                         !Interface quantities
  type(derivative_data_type_df5), dimension(3) :: grad_u, grad_v, grad_w   !Interface gradients of velocities
  type(derivative_data_type_df5), dimension(3) :: grad_rho, grad_p, grad_T !Interface gradients of rho, p, and T
 
  real(p2), dimension(3) :: grad_uL, grad_vL, grad_wL, grad_rL, grad_pL
  real(p2), dimension(3) :: grad_uR, grad_vR, grad_wR, grad_rR, grad_pR
 
  type(derivative_data_type_df5) :: rho, a2   !Interface values for density and (speed of sound)^2
 
  integer  :: ix, iy, iz    !Index to extract the x-, y-, and z-component.
 
  C = C_0   !Parameter for Sutherland's law
  Re_inf = Reynolds   !Free stream Reynolds number
  T_inf = Freestream_Temp   !Free stream temperature
  
     alpha = one
 
        ix = 1
        iy = 2
        iz = 3
 
 ! Left and right states and gradients:
 
        rhoL = ucL(1)
          uL = ucL(2)/ucL(1)
          vL = ucL(3)/ucL(1)
          wL = ucL(4)/ucL(1)
          pL = (gamma-one)*( ucL(5) - half*rhoL*(uL*uL+vL*vL+wL*wL) )
 
     grad_rL = gradwL(:,1)
     grad_uL = gradwL(:,2)
     grad_vL = gradwL(:,3)
     grad_wL = gradwL(:,4)
     grad_pL = gradwL(:,5)
 
        rhoR = ucR(1)
          uR = ucR(2)/ucR(1)
          vR = ucR(3)/ucR(1)
          wR = ucR(4)/ucR(1)
          pR = (gamma-one)*( ucR(5) - half*rhoR*(uR*uR+vR*vR+wR*wR) )
 
     grad_rR = gradwR(:,1)
     grad_uR = gradwR(:,2)
     grad_vR = gradwR(:,3)
     grad_wR = gradwR(:,4)
     grad_pR = gradwR(:,5)
 
 ! Temperature is identical to the speed of sound due to
 ! the nondimensionalization.
 
        TL = gamma*pL/rhoL
        TR = gamma*pR/rhoR
 
 ! Arithmetic averages of velocities and temperature.
 
         u = half*(uL + uR)
         v = half*(vL + vR)
         w = half*(wL + wR)
         T = half*(TL + TR)
 
 ! Sutherland's law in the nondimensional form.
 ! Note: The factor, M_inf/Re_inf, comes from nondimensionalization.
 
        mu = (one+C/T_inf)/(T+C/T_inf)*T**(three*half) * M_inf/Re_inf
 
 ! Damping coefficient, alpha:
 !  (1)alpha=1   gives the central-difference formula in 1D.
 !  (2)alpha=4/3 corresponds to the 4th-order diffusion scheme in 1D.
 
     alpha = four_third
 
 ! Lr = Length scale involving the skewness measure (see Nishikawa AIAA2010-5093)
 ! This is the key quantity for robust and accurate computations on skewed grids.
 
        Lr = half*abs( njk(ix)*ejk(ix) + njk(iy)*ejk(iy) + njk(iz)*ejk(iz) ) * mag_ejk
 
 ! Interface gradients from the derived diffusion scheme (Nishikawa-AIAA2010-5093).
 ! The second term is the damping term.
 
    grad_u = half*( (grad_uR + grad_uL) + alpha/Lr*(uR-uL)*njk )
    grad_v = half*( (grad_vR + grad_vL) + alpha/Lr*(vR-vL)*njk )
    grad_w = half*( (grad_wR + grad_wL) + alpha/Lr*(wR-wL)*njk )
 
 ! The temperature gradient is computed from the interface density and the pressure
 ! gradients. 
 ! Note: T = gamma*p/rho -> grad(T) = gamma*grad(p)/rho - (gamma*p/rho^2)*grad(rho)
 
        rho = half*(rhoR + rhoL)
         a2 = gamma*half*(pR + pL)/rho
 
     grad_rho = half*( (grad_rR + grad_rL) + alpha/Lr*(rhoR-rhoL)*njk )
     grad_p   = half*( (grad_pR + grad_pL) + alpha/Lr*(  pR-pL  )*njk )
 
     grad_T   = ( gamma*grad_p - a2*grad_rho) /rho
 
 ! Interface gradients have been computed: grad(u), grad(v), grad(w), grad(T).
 ! We now evaluate the physical viscous flux with them.
 
 ! Viscous stresses (Stokes' hypothesis is assumed)
 
     tauxx =  mu*(four_third*grad_u(ix) - two_third*grad_v(iy) - two_third*grad_w(iz))
     tauyy =  mu*(four_third*grad_v(iy) - two_third*grad_u(ix) - two_third*grad_w(iz))
     tauzz =  mu*(four_third*grad_w(iz) - two_third*grad_u(ix) - two_third*grad_v(iy))
 
     tauxy =  mu*(grad_u(iy) + grad_v(ix))
     tauxz =  mu*(grad_u(iz) + grad_w(ix))
     tauyz =  mu*(grad_v(iz) + grad_w(iy))
 
     tauyx = tauxy
     tauzx = tauxz
     tauzy = tauyz
 
 ! Heat fluxes: q = - mu*grad(T)/(Prandtl*(gamma-1))
 
        qx = - mu*grad_T(ix)/(Prandtl*(gamma-one))
        qy = - mu*grad_T(iy)/(Prandtl*(gamma-one))
        qz = - mu*grad_T(iz)/(Prandtl*(gamma-one))
 
 ! Normal components
 
     tauxn = tauxx*njk(ix) + tauxy*njk(iy) + tauxz*njk(iz)
     tauyn = tauyx*njk(ix) + tauyy*njk(iy) + tauyz*njk(iz)
     tauzn = tauzx*njk(ix) + tauzy*njk(iy) + tauzz*njk(iz)
        qn =    qx*njk(ix) +    qy*njk(iy) +    qz*njk(iz)
 
 ! Evaluate the viscous flux at the interface
 
     numerical_flux(1) =   zero
     numerical_flux(2) = - tauxn
     numerical_flux(3) = - tauyn
     numerical_flux(4) = - tauzn
     numerical_flux(5) = - (tauxn*u + tauyn*v + tauzn*w) + qn
 
 ! Normal max wave speed
 !    wsn = alpha*(mu/rho*gamma/Prandtl)/Lr
 
  end subroutine viscous_alpha_ddt
 !--------------------------------------------------------------------------------
 
 end module flux_functions_ddt
