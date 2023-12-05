module module_flux
    implicit none
    public :: interface_flux
    public :: perturb_flux

    ! List of available flux functions:
    private :: roe ! <-- Roe flux

    contains
    subroutine interface_flux(u1, u2, gradw1, gradw2, n12, &
                             xc1, yc1, zc1, xc2, yc2, zc2, &
                              xm, ym, zm, phi1, phi2, num_flux, wsn, uR2L, uR2R )
        use module_common_data     , only : p2, zero
        use module_input_parameter , only : inviscid_flux !name of flux specified as input
        use module_input_parameter , only : accuracy_order, navier_stokes
        use module_input_parameter , only : low_mach_correction, eps_weiss_smith, min_ref_vel, pressure_dif_ws
        use module_ccfv_data_soln  , only : u2w, w2u, q2u !variable conversion functions  

        implicit none

        ! Inputs
        real(p2), dimension(5),     intent(in) :: u1, u2            ! Conservative vars (prim for low mach)
        real(p2), dimension(3,5),   intent(in) :: gradw1, gradw2    ! Gradients of primitive vars
        real(p2), dimension(3),     intent(in) :: n12               ! Unit area vector (from c1 to c2)
        real(p2),                   intent(in) :: xc1, yc1, zc1     ! Left cell centroid
        real(p2),                   intent(in) :: xc2, yc2, zc2     ! Right cell centroid
        real(p2),                   intent(in) :: xm, ym, zm        ! Face midpoint
        real(p2),                   intent(in) :: phi1, phi2        ! Limiter

        real(p2), optional,         intent(in) :: uR2L, uR2R        ! preconditioner scaling factor
        ! Output
        real(p2), dimension(5),     intent(out) :: num_flux         ! Output
        real(p2),                   intent(out) :: wsn              ! Max wave speed at face

        ! Local Vars
        real(p2), dimension(5) :: w1, w2 ! Primitive vars at centroids
        real(p2), dimension(5) :: wL, wR, qL, qR ! primitive vars reconstructed to face
        real(p2), dimension(5) :: uL, uR ! conservative vars computed from wL and wR
        real(p2)               :: uR2L_rec, uR2R_rec ! Reconstructed uR2
        real(p2)               :: dp ! difference in pressure, needed for low mach preconditioning
        real(p2), dimension(3) :: ejk
        REAL(P2)               :: mag_ejk ! step difference for VISCOUS U_REF

        if (low_mach_correction) then !(q1 = u1, q2 = u2)
            if (accuracy_order == 2) then
                qL = u1 + phi1 * ( gradw1(1,:)*(xm-xc1) + gradw1(2,:)*(ym-yc1) + gradw1(3,:)*(zm-zc1) ) ! gradq <=> gradw (var) 
                qR = u2 + phi2 * ( gradw2(1,:)*(xm-xc2) + gradw2(2,:)*(ym-yc2) + gradw2(3,:)*(zm-zc2) ) ! u <=> q (vars)
            else
                qL = u1 
                qR = u2
            end if
            uL = q2u(qL)
            uR = q2u(qR)
            ! dp =   abs(qR(1) - qL(1))   !Pressure difference
            ! uR2L_rec = max(sqrt(qL(2)**2 + qL(3)**2 + qL(4)**2),pressure_dif_ws*sqrt(dp/uL(1)), eps_weiss_smith*qL(5), min_ref_vel)
            ! uR2R_rec = max(sqrt(qR(2)**2 + qR(3)**2 + qR(4)**2),pressure_dif_ws*sqrt(dp/uR(1)), eps_weiss_smith*qR(5), min_ref_vel)
            ! if (navier_stokes) then
            !     ! ! uR2R_rec = max(uR2R_rec,)
            !     ! ejk = (/ xc2 - xc1, yc2 - yc1, zc2 - zc1 /)
            !     ! mag_ejk = sqrt(ejk(1)**2 + ejk(2)**2 + ejk(3)**2)

            !     ! uR2l_rec = max(uR2L_rec,mu1/uL(1)/mag_ejk)
            !     ! uR2R_rec = max(uR2R_rec,mu2/uR(2)/mag_ejk)
            ! end if
            ! uR2L_rec = min(uR2L_rec,qL(5))**2 ! cap with the speed of sound
            ! uR2R_rec = min(uR2R_rec,qR(5))**2
            uR2L_rec = uR2L
            uR2R_rec = uR2R
        else
            w1 = u2w(u1)
            w2 = u2w(u2)
            ! Linear reconstruction
            if (accuracy_order == 2) then
                wL = w1 + phi1 * ( gradw1(1,:)*(xm-xc1) + gradw1(2,:)*(ym-yc1) + gradw1(3,:)*(zm-zc1) )
                wR = w2 + phi2 * ( gradw2(1,:)*(xm-xc2) + gradw2(2,:)*(ym-yc2) + gradw2(3,:)*(zm-zc2) )
            else
                wL = w1
                wR = w2
            end if
            ! Flux function uses conservative vars
            uL = w2u(wL)
            uR = w2u(wR)
        end if
        
        !------------------------------------------------------------
        !  (1) Roe flux
        !------------------------------------------------------------
        if(trim(inviscid_flux)=="roe") then
            if (low_mach_correction) then
                    call roe_low_mach(uL,uR,uR2L_rec,uR2R_rec,n12,num_flux,wsn)
            else
                call roe(uL,uR,n12, num_flux,wsn)
            end if
        ! !------------------------------------------------------------
        ! Other fluxes not yet implemneted.
        !------------------------------------------------------------
        else

            write(*,*) " Invalid input for inviscid_flux = ", trim(inviscid_flux)
            write(*,*) " Choose roe or rhll, and try again."
            write(*,*) " ... Stop."
            stop

        endif


    end subroutine interface_flux

    subroutine perturb_flux(u1, u2, gradw1, gradw2, n12, &
                            xc1, yc1, zc1, xc2, yc2, zc2, &
                            xm, ym, zm, phi1, phi2, ivar, perturb_val, side, &
                            num_flux, wsn, uR2L, uR2R  )
        use module_common_data     , only : p2, zero
        use module_input_parameter , only : inviscid_flux !name of flux specified as input
        use module_input_parameter , only : accuracy_order, low_mach_correction
        use module_ccfv_data_soln  , only : u2w, w2u, q2u      !variable conversion functions  

        implicit none

        ! Inputs
        real(p2), dimension(5),     intent(in) :: u1, u2            ! Conservative vars
        real(p2), dimension(3,5),   intent(in) :: gradw1, gradw2    ! Gradients of primitive vars
        real(p2), dimension(3),     intent(in) :: n12               ! Unit area vector (from c1 to c2)
        real(p2),                   intent(in) :: xc1, yc1, zc1     ! Left cell centroid
        real(p2),                   intent(in) :: xc2, yc2, zc2     ! Right cell centroid
        real(p2),                   intent(in) :: xm, ym, zm        ! Face midpoint
        integer,                    intent(in) :: ivar, side        ! Variable to be perturbed
        real(p2),                   intent(in) :: perturb_val       ! Size of perturbation
        real(p2),                   intent(in) :: phi1, phi2        ! Limiter
        real(p2), optional,         intent(in) :: uR2L, uR2R        ! preconditioner scaling factor

        ! Output
        real(p2), dimension(5),     intent(out) :: num_flux         ! Output
        real(p2),                   intent(out) :: wsn              ! Max wave speed at face

        ! Local Vars
        real(p2), dimension(5) :: w1, w2 ! Primitive vars at centroids
        real(p2), dimension(5) :: wL, wR, qL, qR ! primitive vars reconstructed to face
        real(p2), dimension(5) :: uL, uR ! conservative vars computed from wL and wR
        real(p2)               :: uR2L_pert, uR2R_pert
        real(p2)               :: epsilon = 10e-05_p2 ! isn't that just 1e-04?
        

        if (low_mach_correction) then !(q1 = u1, q2 = u2)
            if (accuracy_order == 2) then
                qL = u1 + phi1 * ( gradw1(1,:)*(xm-xc1) + gradw1(2,:)*(ym-yc1) + gradw1(3,:)*(zm-zc1) ) ! gradq <=> gradw (var) 
                qR = u2 + phi2 * ( gradw2(1,:)*(xm-xc2) + gradw2(2,:)*(ym-yc2) + gradw2(3,:)*(zm-zc2) ) ! u <=> q (vars)
            else
                qL = u1 
                qR = u2
            end if
            if (side == 1) then
                qL(ivar) = qL(ivar) + perturb_val
                if (sqrt(qL(2)**2 + qL(3)**2 + qL(4)**2) < epsilon * qL(5)) then
                    uR2L_pert = (epsilon * qL(5))**2
                 elseif (sqrt(qL(2)**2 + qL(3)**2 + qL(4)**2) < qL(5)) then
                    uR2L_pert = qL(2)**2 + qL(3)**2 + qL(4)**2
                 else
                    uR2L_pert = qL(5)**2
                 endif
                uR2R_pert = uR2R
            else if (side == 2) then
                qR(ivar) = qR(ivar) + perturb_val
                if (sqrt(qR(2)**2 + qR(3)**2 + qR(4)**2) < epsilon * qR(5)) then
                    uR2R_pert = (epsilon * qR(5))**2
                 elseif (sqrt(qR(2)**2 + qR(3)**2 + qR(4)**2) < qR(5)) then
                    uR2R_pert = qR(2)**2 + qR(3)**2 + qR(4)**2
                 else
                    uR2R_pert = qR(5)**2
                 endif
                 uR2L_pert = uR2L
            end if
            uL = q2u(qL)
            uR = q2u(qR)
        else
            w1 = u2w(u1)
            w2 = u2w(u2)
            ! Linear reconstruction
            if (accuracy_order == 2) then
                wL = w1 + phi1 * ( gradw1(1,:)*(xm-xc1) + gradw1(2,:)*(ym-yc1) + gradw1(3,:)*(zm-zc1) )
                wR = w2 + phi2 * ( gradw2(1,:)*(xm-xc2) + gradw2(2,:)*(ym-yc2) + gradw2(3,:)*(zm-zc2) )
            else
                wL = w1
                wR = w2
            end if
            ! Flux function uses conservative vars
            uL = w2u(wL)
            uR = w2u(wR)
            if (side == 1) then
                uL(ivar) = uL(ivar) + perturb_val
            else if (side == 2) then
                uR(ivar) = uR(ivar) + perturb_val
            end if
        end if

        

        !------------------------------------------------------------
        !  (1) Roe flux
        !------------------------------------------------------------
        if(trim(inviscid_flux)=="roe") then
            if (low_mach_correction) then
                if (present(uR2L) .and. present(uR2R)) then
                    call roe_low_mach(uL,uR,uR2L_pert,uR2R_pert,n12,num_flux,wsn)
            !         write(*,*) "UR2L = ", uR2L, "UR2R = ", uR2R
            !         write(*,*) "Low Mach:", num_flux(:)
            !         write(*,*) "Low Mach WSN:", wsn
                else
                    write(*,*) "uR2L and/or uR2R not present. Stop!"
                    stop
                end if
            else
                call roe(uL,uR,n12, num_flux,wsn)
            !     write(*,*) " No Prec:", num_flux(:)
            !     write(*,*) " No Mach WSN:", wsn
            !     write(*,*)

            end if

        !------------------------------------------------------------
        ! Other fluxes not yet implemneted.
        !------------------------------------------------------------
        else

            write(*,*) " Invalid input for inviscid_flux = ", trim(inviscid_flux)
            write(*,*) " Choose roe or rhll, and try again."
            write(*,*) " ... Stop."
            stop

        endif


    end subroutine perturb_flux

    subroutine viscous_alpha(ucL,ucR,gradwL,gradwR, njk,ejk,mag_ejk, numerical_flux)

        ! use derivative_data_df5
        use module_input_parameter , only : M_inf, Reynolds, C_0, Freestream_Temp


        implicit none
        integer , parameter :: p2 = selected_real_kind(15)
       
       !Input
        real(p2)                      , dimension(5)  , intent( in) :: ucL     !Left state (conservative)
        real(p2)                      , dimension(5)  , intent( in) :: ucR     !Right state (conservative)
        real(p2)                      , dimension(3,5), intent( in) :: gradwL  !Left gradient (primitive)
        real(p2)                      , dimension(3,5), intent( in) :: gradwR  !Right gradient (primitive)
        real(p2)                      , dimension(3)  , intent( in) :: njk     !Unit directed area vector
        real(p2)                      , dimension(3)  , intent( in) :: ejk     !Unit edge vector
        real(p2)                      ,                 intent( in) :: mag_ejk !Magnitude of the edge vector
       
       !Output
        real(p2)                      , dimension(5), intent(out) :: numerical_flux !Numerical viscous flux
       
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
       
        real(p2)                       ::   uL, uR            ! x-velocity  (Left and Right states)
        real(p2)                       ::   vL, vR            ! y-velocity  (Left and Right states)
        real(p2)                       ::   wL, wR            ! z-velocity  (Left and Right states)
        real(p2)                       :: rhoL, rhoR          ! Density     (Left and Right states)
        real(p2)                       ::   pL, pR            ! Pressure    (Left and Right states)
        real(p2)                       ::   TL, TR            ! Temperature (Left and Right states)
       
        real(p2)                       :: tauxx, tauyy, tauzz !Viscous stresses: diagonal compontens
        real(p2)                       :: tauxy, tauyz, tauzx !Viscous stresses: off-diagonal components
        real(p2)                       :: tauyx, tauzy, tauxz !Viscous stresses: same as above by symmetry
        real(p2)                       :: qx, qy, qz          !Heat flux components
        real(p2)                       :: tauxn, tauyn, tauzn !Normal stresses
        real(p2)                       :: qn                  !Normal heat flux
       
        real(p2)                       :: u, v, w, T, mu                         !Interface quantities
        real(p2)                      , dimension(3) :: grad_u, grad_v, grad_w   !Interface gradients of velocities
        real(p2)                      , dimension(3) :: grad_rho, grad_p, grad_T !Interface gradients of rho, p, and T
       
        real(p2), dimension(3) :: grad_uL, grad_vL, grad_wL, grad_rL, grad_pL
        real(p2), dimension(3) :: grad_uR, grad_vR, grad_wR, grad_rR, grad_pR
       
        real(p2)                       :: rho, a2   !Interface values for density and (speed of sound)^2
       
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
              ! get_viscosity = scaling_factor * ( (one + ( C_0/Freestream_Temp ) )/(T + ( C_0/Freestream_Temp )) ) ** 1.5_p2
              if (isnan(mu)) then 
                write (*,*) "nan value present - press [Enter] to continue"
                read(unit=*,fmt=*)
              end if

            !   mu = ( (one + ( C_0/T_inf ) )/(T + ( C_0/Freestream_Temp )) ) ** 1.5_p2
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
       
        end subroutine viscous_alpha





    subroutine viscous_flux(w1,w2,gradw1,gradw2,T1,T2,gradT1,gradT2,mu1,mu2,n12,xc1,yc1,zc1,xc2,yc2,zc2, visc_flux)
        use module_common_data, only : p2, half, one, zero
        use module_ccfv_data_soln, only : gamma, w2u
        use module_input_parameter, only : Pr


        implicit none

        real(p2), dimension(5),     intent(in) :: w1, w2
        real(p2), dimension(3,5),   intent(in) :: gradw1, gradw2
        real(p2),                   intent(in) :: T1,T2, mu1, mu2
        real(p2), dimension(3),     intent(in) :: gradT1, gradT2
        real(p2), dimension(3),     intent(in) :: n12               ! Unit area vector (from c1 to c2)
        real(p2),                   intent(in) :: xc1, yc1, zc1     ! Left cell centroid
        real(p2),                   intent(in) :: xc2, yc2, zc2     ! Right cell centroid
        real(p2), dimension(5),     INTENT(OUT):: visc_flux
        
        
        ! local vars
        real(p2), dimension(3,3)    :: face_gradw, gradQ, tau
        real(p2), dimension(3)      :: face_gradT, ds, dsds2, delU, w_face
        ! face_gradT = [dTdx dTdy dTdz]
        real(p2)                    :: mu_face, heat_conductivity
        real(p2), dimension(5)      :: u1, u2
        real(p2)                    :: theta_1, theta_2, theta_3

        integer                     :: i, iu
        
        real(p2), dimension(3) :: grad_uL, grad_vL, grad_wL, grad_rL, grad_pL
        real(p2), dimension(3) :: grad_uR, grad_vR, grad_wR, grad_rR, grad_pR

        ! Calculate the face gradients
        ds = (/xc2-xc1, yc2-yc1, zc2-zc1/) ! vector pointing from center of cell 1 to cell 2
        dsds2 = ds/(ds(1)**2 + ds(2)**2 + ds(3)**2) ! ds/ds^2

        grad_rL = gradw1(:,1)
        grad_uL = gradw1(:,2)
        grad_vL = gradw1(:,3)
        grad_wL = gradw1(:,4)
        grad_pL = gradw1(:,5)

        grad_rR = gradw2(:,1)
        grad_uR = gradw2(:,2)
        grad_vR = gradw2(:,3)
        grad_wR = gradw2(:,4)
        grad_pR = gradw2(:,5)
        
        do iu = 1,3
            delU = half * (gradw1(:,iu+1) + gradw2(:,iu+1))
            face_gradw(iu,:) = delU + ((w2(iu+1) - w1(iu+1)) - dot_product(delU,ds) ) * dsds2
            ! delU = delU_bar + [dU - dot(delU_bar,ds)]*ds/ds^2
            ! delU_bar is the arithmetic mean of the left and right gradients.  In order to prevent checkerboarding on certain 
            ! grids the gradient along the vector ds is replaced with a central difference.
        end do
        
        
        delU = half * (gradT1(:) + gradT2(:))
        face_gradT(:) = delU + ((T2 - T1) - dot_product(delU,ds) ) * dsds2
        mu_face = half * (mu1 + mu2) ! mu = M_inf*mu_ND/Re_inf
        tau = compute_tau(face_gradw,mu_face)
        heat_conductivity = mu_face/((gamma - one) * Pr)
        w_face = half * (w1(2:4) + w2(2:4))

        theta_1 = -(w_face(1)*tau(1,1) + w_face(2)*tau(2,1) + w_face(3)*tau(3,1)) + heat_conductivity * face_gradT(1)
        theta_2 = -(w_face(1)*tau(1,2) + w_face(2)*tau(2,2) + w_face(3)*tau(3,2)) + heat_conductivity * face_gradT(2)
        theta_3 = -(w_face(1)*tau(1,3) + w_face(2)*tau(2,3) + w_face(3)*tau(3,3)) + heat_conductivity * face_gradT(3)


        visc_flux(1) = zero
        visc_flux(2) = -(n12(1)*tau(1,1) + n12(2)*tau(1,2) + n12(3)*tau(1,3))
        visc_flux(3) = -(n12(1)*tau(2,1) + n12(2)*tau(2,2) + n12(3)*tau(2,3))
        visc_flux(4) = -(n12(1)*tau(3,1) + n12(2)*tau(3,2) + n12(3)*tau(3,3))
        visc_flux(5) = n12(1)*theta_1  + n12(2)*theta_2  + n12(3)*theta_3

        
    end subroutine viscous_flux

    subroutine viscous_flux_primitive(q1,q2,gradq1,gradq2,mu1,mu2,n12,xc1,yc1,zc1,xc2,yc2,zc2, visc_flux)
        use module_common_data, only : p2, half, one, zero
        use module_ccfv_data_soln, only : gamma, w2u
        use module_input_parameter, only : Pr


        implicit none

        real(p2), dimension(5),     intent(in) :: q1, q2
        real(p2), dimension(3,5),   intent(in) :: gradq1, gradq2
        real(p2),                   intent(in) :: mu1, mu2
        real(p2), dimension(3),     intent(in) :: n12               ! Unit area vector (from c1 to c2)
        real(p2),                   intent(in) :: xc1, yc1, zc1     ! Left cell centroid
        real(p2),                   intent(in) :: xc2, yc2, zc2     ! Right cell centroid
        real(p2), dimension(5),     INTENT(OUT):: visc_flux
        
        
        ! local vars
        real(p2), dimension(3,3)    :: face_gradq, gradQ, tau
        real(p2), dimension(3)      :: face_gradT, ds, dsds2, delU, q_face
        ! face_gradT = [dTdx dTdy dTdz]
        real(p2)                    :: mu_face, heat_conductivity
        real(p2), dimension(5)      :: u1, u2
        real(p2)                    :: theta_1, theta_2, theta_3

        integer                     :: i, iu
        
        real(p2), dimension(3) :: grad_uL, grad_vL, grad_wL, grad_TL, grad_pL
        real(p2), dimension(3) :: grad_uR, grad_vR, grad_wR, grad_TR, grad_pR

        ! Calculate the face gradients
        ds = (/xc2-xc1, yc2-yc1, zc2-zc1/) ! vector pointing from center of cell 1 to cell 2
        dsds2 = ds/(ds(1)**2 + ds(2)**2 + ds(3)**2) ! ds/ds^2

        grad_pL = gradq1(:,1)
        grad_uL = gradq1(:,2)
        grad_vL = gradq1(:,3)
        grad_wL = gradq1(:,4)
        grad_TL = gradq1(:,5)

        grad_pR = gradq2(:,1)
        grad_uR = gradq2(:,2)
        grad_vR = gradq2(:,3)
        grad_wR = gradq2(:,4)
        grad_TR = gradq2(:,5)
        
        do iu = 1,3
            ! delU = delU_bar + [dU - dot(delU_bar,ds)]*ds/ds^2
            ! delU_bar is the arithmetic mean of the left and right gradients.  In order to prevent checkerboarding on certain 
            ! grids the gradient along the vector ds is replaced with a central difference.
            delU = half * (gradq1(:,iu+1) + gradq2(:,iu+1))
            face_gradq(iu,:) = delU + ((q2(iu+1) - q1(iu+1)) - dot_product(delU,ds) ) * dsds2
        end do
        
        
        delU = half * (gradq1(:,5) + gradq2(:,5))
        face_gradT(:) = delU + ((q2(5) - q1(5)) - dot_product(delU,ds) ) * dsds2
        mu_face = half * (mu1 + mu2) ! mu = M_inf*mu_ND/Re_inf
        tau = compute_tau(face_gradq,mu_face)
        heat_conductivity = mu_face/((gamma - one) * Pr)
        q_face = half * (q1(2:4) + q2(2:4)) ! three face velocities

        theta_1 = -(q_face(1)*tau(1,1) + q_face(2)*tau(2,1) + q_face(3)*tau(3,1)) + heat_conductivity * face_gradT(1)
        theta_2 = -(q_face(1)*tau(1,2) + q_face(2)*tau(2,2) + q_face(3)*tau(3,2)) + heat_conductivity * face_gradT(2)
        theta_3 = -(q_face(1)*tau(1,3) + q_face(2)*tau(2,3) + q_face(3)*tau(3,3)) + heat_conductivity * face_gradT(3)


        visc_flux(1) = zero
        visc_flux(2) = -(n12(1)*tau(1,1) + n12(2)*tau(1,2) + n12(3)*tau(1,3))
        visc_flux(3) = -(n12(1)*tau(2,1) + n12(2)*tau(2,2) + n12(3)*tau(2,3))
        visc_flux(4) = -(n12(1)*tau(3,1) + n12(2)*tau(3,2) + n12(3)*tau(3,3))
        visc_flux(5) = n12(1)*theta_1  + n12(2)*theta_2  + n12(3)*theta_3

        
    end subroutine viscous_flux_primitive

    subroutine viscous_bflux(w1,wb,T1,Tb,mu1, xc1,yc1,zc1, xb,yb,zb, bn1,bn2, bface_nrml_mag,unit_face_normal,viscous_flux_boundary)
        use module_common_data, only : p2, half, one, zero, x, y, z
        use module_ccfv_data_soln, only : gamma
        use module_input_parameter, only : Pr
        use module_gewp, only : gewp_solve
        
        implicit none
        
        real(p2), dimension(5),     intent(in) :: w1,wb
        real(p2),                   intent(in) :: T1,Tb,mu1
        real(p2),                   intent(in) :: xc1,yc1,zc1       ! cell 1 center
        real(p2),                   intent(in) :: xb,yb,zb          ! boundary face center
        integer,                    intent(in) :: bn1, bn2
        real(p2),                   intent(in) :: bface_nrml_mag
        real(p2), dimension(3),     intent(in) :: unit_face_normal
        real(p2), dimension(5),     INTENT(OUT):: viscous_flux_boundary

        real(p2), dimension(3)      :: d_Cb, V_parallel, wall_shear, face_gradT
        real(p2)                    :: d_perp, mu_face, heat_conductivity
        real(p2), dimension(3,3)    :: del, del_inv

        integer                     :: os
        

        viscous_flux_boundary(1) = zero
        
        mu_face = mu1 ! just set it to the cell value.  This should be sufficient.  We can be more clever later if we want...
    
        ! Vector from cell center to boundary face center: 
        d_Cb = -(/xc1 - xb, yc1 - yb, zc1 - zb/)
        ! Component of d_Cb perpendicular to boundary face:
        d_perp = dot_product(d_Cb,unit_face_normal)

        ! Calculate the wall shear force:
        wall_shear = -(mu_face*bface_nrml_mag/d_perp) * &
            ((w1(2:4) - wb(2:4)) - (dot_product(w1(2:4) - wb(2:4),unit_face_normal)*unit_face_normal))

        viscous_flux_boundary(2:4) = wall_shear

        ! Calculating wall temperature gradient
                    
        ! Calculate two vectors from boundary face center to boundary nodes.
        del(1,:) = d_Cb
        del(2,:) = (/x(bn1), y(bn1), z(bn1)/) - (/xb,yb,zb/)
        del(3,:) = (/x(bn2), y(bn2), z(bn2)/) - (/xb,yb,zb/)

        call gewp_solve(del, 3, del_inv, os)
        
        face_gradT = matmul(del_inv,(/one - T1, zero, zero/))

        heat_conductivity = mu_face/((gamma - one) * Pr) ! T_wall = T_inf = 1 (T is normalized with T_inf)

        viscous_flux_boundary(5) = heat_conductivity *  dot_product(face_gradT,unit_face_normal)
    end subroutine viscous_bflux

    subroutine compute_viscosity
        use module_common_data, only : p2, one
        use module_ccfv_data_grid, only : ncells
        use module_input_parameter, only : M_inf, Reynolds, C_0, Freestream_Temp, low_mach_correction
        use module_ccfv_data_soln, only : Temp, mu, q
        
        implicit none
        integer :: i
        real(p2) :: scaling_factor

        scaling_factor = M_inf/Reynolds

        do i = 1,ncells
            if (low_mach_correction) then
                mu(i) =scaling_factor*( (one + ( C_0/Freestream_Temp ) )/(q(5,i) + ( C_0/Freestream_Temp )) ) * (q(5,i) ** 1.5_p2)
            else
                mu(i) =scaling_factor*( (one + ( C_0/Freestream_Temp ) )/(Temp(i) + ( C_0/Freestream_Temp )) ) * (Temp(i) ** 1.5_p2)
            endif
        end do
        
    end subroutine compute_viscosity

    real(p2)function get_viscosity(T)
        use module_common_data, only : p2, one
        use module_ccfv_data_grid, only : ncells
        use module_input_parameter, only : M_inf, Reynolds, C_0, Freestream_Temp
        use module_ccfv_data_soln, only : Temp, mu
        
        
        implicit none
        real(p2), intent(in) :: T ! Temperature
        integer :: i
        real(p2) :: scaling_factor

        scaling_factor = M_inf/Reynolds

        get_viscosity = scaling_factor * ( (one + ( C_0/Freestream_Temp ) )/(T + ( C_0/Freestream_Temp )) ) ** 1.5_p2
        
        
    end function get_viscosity

    function compute_tau(gradU_face,mu_face) ! lol this is wrong...
          use module_common_data, only : p2, third, two
          use module_ccfv_data_grid, only : ncells

          implicit none
          real(p2), dimension(3,3), intent(in) :: gradU_face
          real(p2), intent(in)                 :: mu_face
          real(p2), dimension(3,3)             :: compute_tau
          
          integer :: i
          real(p2) :: divergence

          
          
          divergence = gradU_face(1,1) + gradU_face(2,2) + gradU_face(3,3)
          compute_tau(1,1) = two*mu_face*(gradU_face(1,1) - third * (divergence)) ! Tau_xx
          compute_tau(2,2) = two*mu_face*(gradU_face(2,2) - third * (divergence)) ! Tau_yy
          compute_tau(3,3) = two*mu_face*(gradU_face(3,3) - third * (divergence)) ! Tau_zz
          compute_tau(1,2) = mu_face * (gradU_face(2,1) + gradU_face(1,2))           ! Tau_xy = mu(du/dy + dv/dx)
          compute_tau(2,1) = compute_tau(1,2)
          compute_tau(1,3) = mu_face * (gradU_face(3,1) + gradU_face(1,3))           ! Tau_xz = mu(du/dz + dw/dx)
          compute_tau(3,1) = compute_tau(1,3)
          compute_tau(2,3) = mu_face * (gradU_face(3,2) + gradU_face(2,3))           ! Tau_xy = mu(du/dy + dv/dx)
          compute_tau(3,2) = compute_tau(2,3)
          
    end Function compute_tau

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
    !*                    wsn = maximum wave speed (eigenvalue)
    !*
    !* ------------------------------------------------------------------------------
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
    ! 
    !
    !
    ! Note: This is currently unaltered from the edu_euler roe function (other than some formatting changes).  
    ! I will at some point optimie this.  But right now I'm just looking to get a working code...
    subroutine roe(ucL, ucR, njk, num_flux,wsn)

        use module_ccfv_data_soln , only : gamma
        use module_input_parameter, only : eig_limiting_factor
       
        implicit none
       
        integer , parameter :: p2 = selected_real_kind(15) ! Double precision
       
        !Input
        real(p2), dimension(5), intent( in) :: ucL !Left  state in conservative variables.
        real(p2), dimension(5), intent( in) :: ucR !Right state in conservative variables.
        real(p2), dimension(3), intent( in) :: njk
       
        !Output
        real(p2), dimension(5)  , intent(out) :: num_flux !Numerical viscous flux
        real(p2),                 intent(out) :: wsn      !Max wave speed
       
        !Some constants
        real(p2) ::  zero = 0.0_p2
        real(p2) ::   one = 1.0_p2
        real(p2) ::   two = 2.0_p2
        real(p2) ::  half = 0.5_p2
       
        !Local variables
        !            L = Left
        !            R = Right
        ! No subscript = Roe average
       
        real(p2) :: nx, ny, nz             ! Normal vector components
        real(p2) :: uL, uR, vL, vR, wL, wR ! Velocity components.
        real(p2) :: rhoL, rhoR, pL, pR     ! Primitive variables.
        real(p2) :: qnL, qnR               ! Normal velocities
        real(p2) :: aL, aR, HL, HR         ! Speed of sound, Total enthalpy
        real(p2), dimension(5)   :: fL     ! Physical flux evaluated at ucL
        real(p2), dimension(5)   :: fR     ! Physical flux evaluated at ucR
       
        real(p2) :: RT                     ! RT = sqrt(rhoR/rhoL)
        real(p2) :: rho,u,v,w,H,a,qn       ! Roe-averages
       
        real(p2) :: drho, dqn, dp          ! Differences in rho, qn, p, e.g., dp=pR-pL
        real(p2), dimension(4) :: LdU      ! Wave strengths = L*(UR-UL)
        real(p2) :: du, dv, dw             ! Velocity differences
        real(p2), dimension(4)   :: ws     ! Wave speeds
        real(p2), dimension(4)   :: dws    ! Width of a parabolic fit for entropy fix
        real(p2), dimension(5,4) :: R      ! Right-eigenvector matrix
        real(p2), dimension(5)   :: diss   ! Dissipation term
       
        integer  :: i
       
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
             aL = sqrt(gamma*pL/rhoL)
             HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL+wL*wL)
       
        !  Right state
       
           rhoR = ucR(1)
             uR = ucR(2)/ucR(1)
             vR = ucR(3)/ucR(1)
             wR = ucR(4)/ucR(1)
            qnR = uR*nx + vR*ny + wR*nz
             pR = (gamma-one)*( ucR(5) - half*rhoR*(uR*uR+vR*vR+wR*wR) )
             aR = sqrt(gamma*pR/rhoR)
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
        
           RT = sqrt(rhoR/rhoL)
          rho = RT*rhoL                                        !Roe-averaged density
            u = (uL + RT*uR)/(one + RT)                        !Roe-averaged x-velocity
            v = (vL + RT*vR)/(one + RT)                        !Roe-averaged y-velocity
            w = (wL + RT*wR)/(one + RT)                        !Roe-averaged z-velocity
            H = (HL + RT*HR)/(one + RT)                        !Roe-averaged total enthalpy
            a = sqrt( (gamma-one)*(H-half*(u*u + v*v + w*w)) ) !Roe-averaged speed of sound
           qn = u*nx + v*ny + w*nz                             !Roe-averaged face-normal velocity
       
        !Wave Strengths
       
          drho = rhoR - rhoL !Density difference
            dp =   pR - pL   !Pressure difference
           dqn =  qnR - qnL  !Normal velocity difference
       
        LdU(1) = (dp - rho*a*dqn )/(two*a*a) !Left-moving acoustic wave strength
        LdU(2) = (dp + rho*a*dqn )/(two*a*a) !Right-moving acoustic wave strength
        LdU(3) =  drho - dp/(a*a)            !Entropy wave strength
        LdU(4) = rho                         !Shear wave strength (not really, just a factor)
       
        !Absolute values of the wave Speeds
    
        ws(1) = abs(qn-a) !Left-moving acoustic wave
        ws(2) = abs(qn+a) !Right-moving acoustic wave
        ws(3) = abs(qn)   !Entropy wave
        ws(4) = abs(qn)   !Shear waves
       
        ! Harten's Entropy Fix JCP(1983), 49, pp357-393. This is typically applied
        ! only for the nonlinear fields (k=1 and 3), but here it is applied to all
        ! for robustness, avoiding vanishing wave speeds by making a parabolic fit
        ! near ws = 0 for all waves.
        ! 02-27-2018: The limiting can be too much for the shear wave and entropy wave.
        !             Flat plate calculation shows that applying it to all contaminates
        !             the solution significantly. So, apply only to the nonlinear waves,
        !             or apply very small limiting to entropy and shear waves.
        !
        ! Note: ws(1) and ws(2) are the nonlinear waves.
       
        do i = 1, 4
            dws(i) = eig_limiting_factor(i)*a
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
       
        ! Right-moving acoustic wave
        R(1,2) = one
        R(2,2) = u + a*nx
        R(3,2) = v + a*ny
        R(4,2) = w + a*nz
        R(5,2) = H + a*qn
       
        ! Entropy wave
        R(1,3) = one
        R(2,3) = u
        R(3,3) = v 
        R(4,3) = w
        R(5,3) = half*(u*u + v*v + w*w)
       
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
        ! write(*,*) "Normal:", diss
        num_flux = half * (fL + fR - diss)
       
        ! Max wave speed normal to the face:
                    wsn = abs(qn) + a
       
        end subroutine roe
        !--------------------------------------------------------------------------------

        subroutine roe_low_mach(ucL, ucR, uR2L, uR2R, njk, num_flux,wsn)
            ! Modified form of Roe's difference splitting method based off of https://doi.org/10.2514/3.12946

            use module_common_data, only : p2, one, two, half, zero
            use module_ccfv_data_soln , only : gamma
            use module_input_parameter, only : eig_limiting_factor, solver_type, entropy_fix

            implicit none

            !Input
            real(p2), dimension(5), intent( in) :: ucL !Left  state in conservative variables.
            real(p2), dimension(5), intent( in) :: ucR !Right state in conservative variables.
            real(p2),               intent( in) :: uR2L, uR2R ! Left and right scaling terms
            real(p2), dimension(3), intent( in) :: njk

            !Output
            real(p2), dimension(5)  , intent(out) :: num_flux !Numerical viscous flux
            real(p2),                 intent(out) :: wsn      !Max wave speed

            !Local variables
            !            L = Left
            !            R = Right
            ! No subscript = Roe average

            real(p2) :: nx, ny, nz             ! Normal vector components
            real(p2) :: uL, uR, vL, vR, wL, wR ! Velocity components.
            real(p2) :: rhoL, rhoR, pL, pR     ! Primitive variables.
            real(p2) :: qnL, qnR               ! Normal velocities
            real(p2) :: aL, aR, HL, HR         ! Speed of sound, Total enthalpy
            real(p2), dimension(5)   :: fL     ! Physical flux evaluated at ucL
            real(p2), dimension(5)   :: fR     ! Physical flux evaluated at ucR

            real(p2) :: RT                     ! RT = sqrt(rhoR/rhoL)
            real(p2) :: rho,u,v,w,H,a,qn,p     ! Roe-averages

            real(p2) :: drho, dqn, dp          ! Differences in rho, qn, p, e.g., dp=pR-pL
            real(p2) :: drhou, drhov, drhow, drhoE ! Differences in conserved vars
            real(p2) :: absU, ws2, ws3         ! accoustic wave speeds |u'+c'| and |u'-c'|
            real(p2) :: uprime, cprime, alpha, beta, uR2
            real(p2) :: cstar, Mstar, delu, delp
            real(p2), dimension(4) :: LdU      ! Wave strengths = L*(UR-UL)
            real(p2) :: du, dv, dw             ! Velocity differences
            real(p2), dimension(4)   :: ws     ! Wave speeds
            real(p2), dimension(4)   :: dws    ! Width of a parabolic fit for entropy fix
            real(p2), dimension(5,4) :: R      ! Right-eigenvector matrix
            real(p2), dimension(5)   :: diss   ! Dissipation term
           
            real(p2), dimension(3,5) :: diss_cartesian

            integer  :: i
           
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
                 aL = sqrt(gamma*pL/rhoL)
                 HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL+wL*wL)
           
            !  Right state
           
               rhoR = ucR(1)
                 uR = ucR(2)/ucR(1)
                 vR = ucR(3)/ucR(1)
                 wR = ucR(4)/ucR(1)
                qnR = uR*nx + vR*ny + wR*nz
                 pR = (gamma-one)*( ucR(5) - half*rhoR*(uR*uR+vR*vR+wR*wR) )
                 aR = sqrt(gamma*pR/rhoR)
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
           
            !First compute the arithmetic-averaged quantities
            
              rho = half * (rhoR + rhoL)                           !Arithemtic-averaged density
                u = half * (uL   + uR  )                           !Arithemtic-averaged x-velocity
                v = half * (vL   + vR  )                           !Arithemtic-averaged y-velocity
                w = half * (wL   + wR  )                           !Arithemtic-averaged z-velocity
                H = half * (HL   + HR  )                           !Arithemtic-averaged total enthalpy
                p = half * (pL   + pR  )                           !Arithmetic-averaged pressure
                a = sqrt( (gamma-one)*(H-half*(u*u + v*v + w*w)) ) !Arithemtic-averaged speed of sound
               qn = u*nx + v*ny + w*nz                             !Arithemtic-averaged face-normal velocity
              uR2 = half * (uR2L + uR2R)                           !Arithmetic-averaged scaling term


            ! Roe averages
            ! RT = sqrt(rhoR/rhoL)
            ! rho = RT*rhoL                                        !Roe-averaged density
            ! u = (uL + RT*uR)/(one + RT)                        !Roe-averaged x-velocity
            ! v = (vL + RT*vR)/(one + RT)                        !Roe-averaged y-velocity
            ! w = (wL + RT*wR)/(one + RT)                        !Roe-averaged z-velocity
            ! H = (HL + RT*HR)/(one + RT)                        !Roe-averaged total enthalpy
            ! a = sqrt( (gamma-one)*(H-half*(u*u + v*v + w*w)) ) !Roe-averaged speed of sound
            ! qn = u*nx + v*ny + w*nz                             !Roe-averaged face-normal velocity
            ! uR2 = a**2

            !Wave Strengths
           
              drho = rhoR - rhoL !Density difference
                dp =   pR - pL   !Pressure difference
               dqn =  qnR - qnL  !Normal velocity difference (= delta_u in Weiss and Smith)

               drhou = ucR(2) - ucL(2)
               drhov = ucR(3) - ucL(3)
               drhow = ucR(4) - ucL(4)
               drhoE = ucR(5) - ucL(5)
               absU  = abs(qn) ! wave speed one
            !    ! Entropy fix? I guess not....... Need to investigate more
            !     if (absU < 0.1_p2 * a) then
            !         ! if ( ws(i) < dws(i) ) ws(i) = half * ( ws(i)*ws(i)/dws(i)+dws(i) )
            !         absU  = half * ( (absU**2)/(0.2_p2*a) + (0.2_p2*a) ) 
            !     end if

                beta = rho/(gamma*p)
               alpha = half * (one-beta*uR2)
            !    alpha = zero
              cprime = sqrt((alpha**2) * (qn**2) + uR2)
              uprime = qn * (one - alpha)

              ws2 = abs(uprime + cprime)
              ws3 = abs(uprime - cprime)

            !   if (solver_type == 'implicit') then
                ! Entropy fix, only for the implicit solver...
            if (entropy_fix == 'harten') then
                dws(:) = eig_limiting_factor(1:4)*a
                if (absU < dws(3))  absU = half * ( (absU**2)/(dws(3)*a) + (dws(3)*a) )
                if (ws2  < dws(1))  ws2  = half * ( (ws2**2) /(dws(1)*a) + (dws(1)*a) )
                if (ws3  < dws(2))  ws3  = half * ( (ws3**2) /(dws(2)*a) + (dws(2)*a) )
            elseif (entropy_fix == 'mavriplis') then
                ! https://doi.org/10.2514/6.2007-3955
                absU = max(absU,eig_limiting_factor(3)*ws2)
                ! ws2 = max(ws2,eig_limiting_factor(1)*ws2) ! not needed
                ws3 = max(ws3,eig_limiting_factor(1)*ws2)
            endif

              cstar = half * (ws2 + ws3)
              Mstar = half * (ws2 - ws3) / cprime


                delu = Mstar*dqn + (cstar - (one-two*alpha)*absU - alpha*qn*Mstar)*(dp/(rho*uR2))
                delp = Mstar*dp  + (cstar - absU + alpha*qn*Mstar) * rho * dqn

            diss_cartesian(:,1) = (absU * drho + delu * rho) * njk
            diss_cartesian(:,2) = (absU * drhou + delu * rho * u) * njk
            diss_cartesian(:,3) = (absU * drhov + delu * rho * v) * njk
            diss_cartesian(:,4) = (absU * drhow + delu * rho * w) * njk
            diss_cartesian(:,5) = (absU * drhoE + delu * rho * H) * njk

            diss = matmul(transpose(diss_cartesian),njk)  ! compute the component normal to the face to get the actual flux

            ! now add the final term
            ! diss(1) += 0
            diss(2) = diss(2) + delp * nx
            diss(3) = diss(3) + delp * ny
            diss(4) = diss(4) + delp * nz
            diss(5) = diss(5) + delp * qn

            ! This is the numerical flux: Roe flux = 1/2 *[  Fn(UL)+Fn(UR) - GAMMA|An|(WR-WL) ]
            ! write(*,*) "Low-Ma:", diss
            num_flux = half * (fL + fR - diss)
           
            ! Max wave speed normal to the face:
            ! wsn = max(abs(uprime) + cprime,abs(qn)+a)
            wsn = abs(uprime) + cprime
            end subroutine roe_low_mach
end module module_flux