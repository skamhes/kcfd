module module_flux_jac_interface

    implicit none
    public :: interface_jac     !  compute a flux Jacobian at interior face.
    public :: interface_jac_visc ! compute flux jacobian at interior face for viscous flux.
contains

    subroutine interface_jac(wj, wk, njk, dFnduL, dFnduR, uR2L, uR2R )
        use derivative_data_df5
        use module_common_data, only : p2, zero, one, half
        use module_input_parameter, only : inviscid_jac, low_mach_correction
        use flux_functions_ddt    , only :      roe_ddt, &
                                            rusanov_ddt, &
                                                hll_ddt, &
                                               rhll_ddt, &
                                       roe_low_mach_ddt, &
                                       roe_prim_low_mach_ddt
        use module_ccfv_data_soln , only : w2u, u2q

        implicit none

        real(p2), dimension(5), intent(in) :: wj, wk  ! w from cell(j) and neighbor (k)
        real(p2), dimension(3), intent(in) :: njk

        real(p2), dimension(5,5), intent(out) :: dFnduL, dFnduR

        ! Local vavrs
        real(p2), dimension(5)      :: wL, wR, uL, uR
        real(p2), dimension(5,5)    :: dfndu
        real(p2), dimension(5)      :: dummy5
        real(p2)                    :: wsn
        
        real(p2), optional,         intent(in) :: uR2L, uR2R        ! preconditioner scaling factor
        integer :: i, ii
        type(derivative_data_type_df5), dimension(5) :: uL_ddt, uR_ddt, qL_ddt, qR_ddt
        type(derivative_data_type_df5)               :: uR2L_ddt, uR2R_ddt

        ! ! debugging
        ! real(p2) :: H,rho_T,rho_p, theta
        ! real(p2) :: gamma = 1.4_p2
        ! real(p2), dimension(5,5) :: dudq, dfdq

        jac_L_R : do i = 1,2
            if (low_mach_correction) then
                qL_ddt = wj
                qR_ddt = wk
                if (i == 1) then
                    ! Using derivative for uL_ddt
                    call ddt_seed(qL_ddt)
                    uL_ddt = q2u_ddt(qL_ddt)
                    uR_ddt = q2u_ddt(qR_ddt)
                else ! i = 2
                    ! Using derivative for uR_ddt
                    call ddt_seed(qR_ddt)
                    uL_ddt = q2u_ddt(qL_ddt)
                    uR_ddt = q2u_ddt(qR_ddt)
                end if

            else
                ! No reconstruction
                wL = wj
                wR = wk
                ! convert to conservative variables in ddt
                uL_ddt = w2u(wL)
                uR_ddt = w2u(wR)

                ! determine which state the flux derivative is computed wrt
                if (i == 1) then
                    ! Using derivative for uL_ddt
                    call ddt_seed(uL_ddt)
                else ! i = 2
                    ! Using derivative for uR_ddt
                    call ddt_seed(uR_ddt)
                end if
            end if

            !---------------------------------------------------------------------------------
            !---------------------------------------------------------------------------------
            !  Compute inviscid Jacobian by 3D flux subroutines.
            !
            !  Note: These flux subroutines are written based on automatic differentiation,
            !        and thus return the flux and its derivative. Here, we only want the
            !        derivative. So, the flux is received in 'dummy5', and not used.
            !
            !  Note: Input variables to each flux function must be derivative-data-type (ddt),
            !        which contains derivative information as defined in the module
            !        'derivative_data_df5'.
            !
            !---------------------------------------------------------------------------------
            !---------------------------------------------------------------------------------

            !------------------------------------------------------------
            !  (1) Roe flux
            !------------------------------------------------------------
            if(trim(inviscid_jac)=="roe") then
                if (low_mach_correction) then
                    call roe_low_mach_ddt(uL_ddt,uR_ddt,uR2L,uR2R,njk,dummy5,dfndu,wsn)
                else
                    call roe_ddt(uL_ddt,uR_ddt,njk, dummy5,dfndu,wsn)
                end if
            !------------------------------------------------------------
            !  (2) Rusanov flux
            !------------------------------------------------------------
            elseif(trim(inviscid_jac)=="rusanov") then
                call rusanov_ddt(uL_ddt,uR_ddt,njk, dummy5,dfndu,wsn)
            !------------------------------------------------------------
            !  (3) HLL flux
            !------------------------------------------------------------
            elseif(trim(inviscid_jac)=="hll") then
                call hll_ddt(uL_ddt,uR_ddt,njk, dummy5,dfndu,wsn)
            !------------------------------------------------------------
            !  (4) RHLL flux: the last argumant -> exact_jac = .false.
            !                  so that the jac = a1*HLL_jac+a2*Roe_jac
            !                   with a1 and a2 not differentiated.
            !------------------------------------------------------------
            elseif(trim(inviscid_jac)=="rhll") then
                call rhll_ddt(uL_ddt,uR_ddt,njk, dummy5,dfndu,wsn,.false.)
            !------------------------------------------------------------
            !  Others...
            !------------------------------------------------------------
            else
                write(*,*) " Invalid input for inviscid_jac = ", trim(inviscid_jac)
                write(*,*) " Choose roe or rhll, and try again."
                write(*,*) " ... Stop."
                stop
            endif
            if (i==1) then
                dFnduL = dfndu
            else
                dFnduR = dfndu
            endif
        end do jac_L_R
        
    end subroutine interface_jac

    subroutine interface_viscous_alpha_jacobian(uL, uR, gradwL, gradwR, &
        njk, ejk, mag_ejk, dFnduL, dFnduR )

        use module_common_data , only : p2
        use flux_functions_ddt    , only : viscous_alpha_ddt
        use derivative_data_df5


        implicit none
        real(p2)                      , dimension(5)  , intent(in) :: uL, uR
        real(p2)                      , dimension(5,3), intent(in) :: gradwL, gradwR
        real(p2)                      , dimension(3)  , intent(in) :: njk
        real(p2)                      , dimension(3)  , intent(in) :: ejk
        real(p2)                      ,                 intent(in) :: mag_ejk

        real(p2)                      , dimension(5,5), INTENT(OUT):: dFnduL,dFnduR

        integer :: i, ii, jj
        type(derivative_data_type_df5), dimension(5) :: uL_ddt, uR_ddt, flux

        jac_L_R : do i = 1,2
            ! No reconstruction
            ! wL = wj
            ! wR = wk
            ! convert to conservative variables in ddt
            uL_ddt = uL
            uR_ddt = uR

            ! determine which state the flux derivative is computed wrt
            if (i == 1) then
                ! Using derivative for uL_ddt
                call ddt_seed(uL_ddt)
            else ! i = 2
                ! Using derivative for uR_ddt
                call ddt_seed(uR_ddt)
            end if

            call viscous_alpha_ddt(uL_ddt, uR_ddt, gradwL, gradwR, njk,ejk,mag_ejk,  flux)

            if (i==1) then
                do ii = 1,5
                    do jj = 1,5
                        dFnduL(ii,jj) = flux(i)%df(jj)
                    end do
                end do
            else
                do ii = 1,5
                    do jj = 1,5
                        dFnduR(ii,jj) = flux(i)%df(jj)
                    end do
                end do
            endif
        end do jac_L_R
    end subroutine interface_viscous_alpha_jacobian

    subroutine interface_jac_visc(wL, wR, TL, TR, muL, muR, gradVelL, gradVelR, njk, delta_s, dFnduL, dFnduR)

        ! This subroutine builds the interface flux jacobian as outlined in section 4.12 of I do like CFD.
        ! note: only spatial gradients normal to the face are evaluated any gradients wrt global coordinates
        ! (i.e. d/dx, d/dy, and d/dz) are considered to be zero.  While this is not strictly true this is acceptable as
        ! the jacobian only drives convergence (helps pick out the next step).  Convergence will be slower but that's the tradeoff

        
        use module_common_data, only : p2, one, half, two, three, zero
        use module_input_parameter, only : inviscid_jac
        use module_input_parameter, only : M_inf, Reynolds, C_0, Freestream_Temp, Pr
        use module_ccfv_data_soln , only : w2u, gamma

        implicit none

        real(p2), dimension(5), intent(in) :: wL, wR  ! w from cell(j) and neighbor (k)
        real(p2),               intent(in) :: TL,TR,muL,muR ! left and right values for Temperature and vicosity
        real(p2),dimension(3,3),intent(in) :: gradVelL, gradVelR ! gradients of the primative velocity terms
        real(p2), dimension(3), intent(in) :: delta_s ! Vector pointing between cell centers
        real(p2), dimension(3), intent(in) :: njk

        real(p2), dimension(5,5), intent(out) :: dFnduL, dFnduR

        ! Local vavrs
        real(p2), dimension(5,5)    :: dFndw, dwdu
        real(p2), dimension(5)      :: dummy5
        real(p2)                    :: mag_delta_s

        real(p2), dimension(3,3)    :: s_bar
        real(p2), dimension(3)      :: s_bar_n
        real(p2)                    :: s_bar_nv

        real(p2)                    :: Temp_face, mu_face
        real(p2), dimension(3,3)    :: grad_face
        real(p2), dimension(3)      :: dsds2, delU, w_face
        real(p2)                    :: dmudT_bar, sign, scaling_factor, dT_bardrho,dT_bardp
        real(p2)                    :: rho, pressure ! Left or right states of primative vars
        real(p2), dimension(3)      :: vel
        integer :: i, iu
        

        ! Calculate velocity gradients at the face
        dsds2 = delta_s/(delta_s(1)**2 + delta_s(2)**2 + delta_s(3)**2) ! ds/ds^2
        do iu = 1,3
            delU = half * (gradVelL(:,iu) + gradVelR(:,iu))
            grad_face(iu,:) = delU + ((wR(iu+1) - wL(iu+1)) - dot_product(delU,delta_s) ) * dsds2
            ! delU = delU_bar + [dU - dot(delU_bar,ds)]*ds/ds^2
            ! delU_bar is the arithmetic mean of the left and right gradients.  In order to prevent checkerboarding on certain 
            ! grids the gradient along the vector ds is replaced with a central difference.
        end do

        ! Now do the same with temperature
        mu_face = half * (muL + muR) ! mu = M_inf*mu_ND/Re_inf
        Temp_face = half * (TL + TR) 
        w_face = half * (wL(2:4) + wR(2:4))
        
        ! Build the s_bar arrays
        ! diagonal terms
        s_bar(1,1) = (two/three) * (2*grad_face(1,1) - grad_face(2,2) - grad_face(3,3))
        s_bar(2,2) = (two/three) * (2*grad_face(2,2) - grad_face(1,1) - grad_face(3,3))
        s_bar(3,3) = (two/three) * (2*grad_face(3,3) - grad_face(2,2) - grad_face(1,1))
        ! off-diagonal terms (are symmetric)
        s_bar(1,2) = grad_face(1,2) + grad_face(2,1)
        s_bar(2,1) = s_bar(1,2)
        s_bar(1,3) = grad_face(1,3) + grad_face(3,1)
        s_bar(3,1) = s_bar(1,3)
        s_bar(2,3) = grad_face(2,3) + grad_face(3,2)
        s_bar(3,2) = s_bar(2,3)

        ! Now build the s_bar_n vector
        s_bar_n = matmul(s_bar, njk)

        ! And last the s_bar_nv scalar
        s_bar_nv = dot_product(s_bar_n, w_face)

        mag_delta_s = sqrt(delta_s(1)**2 + delta_s(2)**2 + delta_s(3)**2)
        ! Now we can build the jacobian arrays.

        scaling_factor = M_inf/Reynolds

        do i = 1,2
            ! i = 1 => Left
            ! i = 2 => Right
            
            ! First initialize:
            dFndw = zero
            
            ! Calculate some regualr values
            dmudT_bar = (half * (mu_face/(Temp_face + C_0/Freestream_Temp)) * &
                                          (one + three * (C_0/Freestream_Temp)/Temp_face))

            if (i ==1 ) then ! Left
                rho = wL(1)
                vel = wL(2:4)
                pressure = wL(5)
                sign = -one
            else ! Right
                rho = wR(1)
                vel = wR(2:4)
                pressure = wR(5)
                sign = one
            end if

            dT_bardrho = -pressure * gamma / (two * rho**2)
            dT_bardp   = half * gamma / rho

            ! Row 1 = 0

            ! Row 2
            dFndw(2,1) = -dmudT_bar * dT_bardrho * s_bar_n(1)
            dFndw(2,2) = -mu_face * (sign/mag_delta_s)
            ! dFndw(2,3) = zero ! no normal derivative
            ! dFndw(2,4) = zero ! no normal derivative
            dFndw(2,5) = -dmudT_bar * dT_bardp * s_bar_n(1)

            ! Row 3
            dFndw(3,1) = -dmudT_bar * dT_bardrho * s_bar_n(2)
            ! dFndw(3,2) = zero ! no normal derivative
            dFndw(3,3) = -mu_face * (sign/mag_delta_s)
            ! dFndw(3,4) = zero ! no normal derivative
            dFndw(3,5) = -dmudT_bar * dT_bardp * s_bar_n(2)

            ! Row 4
            dFndw(4,1) = -dmudT_bar * dT_bardrho * s_bar_n(3)
            ! dFndw(4,2) = zero ! no normal derivative
            ! dFndw(4,3) = zero ! no normal derivative
            dFndw(4,4) = -mu_face * (sign/mag_delta_s)
            dFndw(4,5) = -dmudT_bar * dT_bardp * s_bar_n(3)

            ! Row 5 
            dFndw(5,1) = -dmudT_bar * dT_bardrho * s_bar_nv & ! dt_n/dW
                -(one/(Pr*(gamma-one))) * (dmudT_bar*dT_bardrho*(TR-TL)+ mu_face * two * sign * dT_bardrho)/mag_delta_s ! dq_n/dW
            dFndw(5,2) = -mu_face * (half * s_bar_n(1) + w_face(1) * (sign/mag_delta_s))
            dFndw(5,3) = -mu_face * (half * s_bar_n(2) + w_face(2) * (sign/mag_delta_s))
            dFndw(5,4) = -mu_face * (half * s_bar_n(3) + w_face(3) * (sign/mag_delta_s))
            dFndw(5,5) = -dmudT_bar * dT_bardp * s_bar_nv & ! dt_n/dW
                -(one/(Pr*(gamma-one))) * (dmudT_bar*dT_bardp*(TR-TL) + mu_face * two * sign * dT_bardp)/mag_delta_s ! dq_n/dW
            
            ! Build dW/dU
            dwdu = zero

            ! Row 1
            dwdu(1,1) = one

            ! Row 2
            dwdu(2,1) = -vel(1)/rho
            dwdu(2,2) = one/rho

            ! Row 3
            dwdu(3,1) = -vel(2)/rho
            dwdu(3,3) = one/rho
            
            ! Row 4
            dwdu(4,1) = -vel(3)/rho
            dwdu(4,4) = one/rho

            ! Row 5
            dwdu(5,1) = half*(gamma - one) * (vel(1)**2 + vel(2)**2 + vel(3)**2)
            dwdu(5,2) = -(gamma-1)*vel(1)
            dwdu(5,3) = -(gamma-1)*vel(2)
            dwdu(5,4) = -(gamma-1)*vel(3)
            dwdu(5,5) =  (gamma-1)

            if (i == 1) then
                dFnduL = matmul(dFndw,dwdu)
            else
                dFnduR = matmul(dFndw,dwdu)
            end if
            
        end do

    end subroutine interface_jac_visc

    !********************************************************************************
    ! Compute U from W (ddt version)
    !
    ! ------------------------------------------------------------------------------
    !  Input:  q =    primitive variables (  p,     u,     v,     w,     T)
    ! Output:  u = conservative variables (rho, rho*u, rho*v, rho*w, rho*E)
    ! ------------------------------------------------------------------------------
    !
    ! Note: rho*E = p/(gamma-1) + rho*0.5*(u^2 + v^2 + w^2)
    !       rho   = p/(gamma*T)
    !********************************************************************************
    function q2u_ddt(q_in) result(u_out)

        use module_common_data, only : p2, one, half
        use module_ccfv_data_soln, only : gamma
        use derivative_data_df5

        implicit none
    
        type(derivative_data_type_df5), dimension(5), intent(in) :: q_in ! input
        type(derivative_data_type_df5), dimension(5)             :: u_out !output
    
        u_out(1) = q_in(1)*gamma / q_in(5)
        u_out(2) = u_out(1)*q_in(2)
        u_out(3) = u_out(1)*q_in(3)
        u_out(4) = u_out(1)*q_in(4)
        u_out(5) = q_in(1)/(gamma-one)+half*u_out(1)*(q_in(2)*q_in(2)+q_in(3)*q_in(3)+q_in(4)*q_in(4))
    
    end function q2u_ddt
end module module_flux_jac_interface