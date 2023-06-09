module module_flux_jac_interface

    implicit none
    public :: interface_jac     !  compute a flux Jacobian at interior face.
contains

    subroutine interface_jac(wj, wk, njk, dFnduL, dFnduR)
        use derivative_data_df5
        use module_common_data, only : p2
        use module_input_parameter, only : inviscid_jac
        use flux_functions_ddt    , only :      roe_ddt, &
                                            rusanov_ddt, &
                                                hll_ddt, &
                                               rhll_ddt
        use module_ccfv_data_soln , only : w2u

        implicit none

        real(p2), dimension(5), intent(in) :: wj, wk  ! w from cell(j) and neighbor (k)
        real(p2), dimension(3), intent(in) :: njk

        real(p2), dimension(5,5), intent(out) :: dFnduL, dFnduR

        ! Local vavrs
        real(p2), dimension(5)      :: wL, wR
        real(p2), dimension(5,5)    :: dfndu
        real(p2), dimension(5)      :: dummy5
        real(p2)                    :: wsn

        integer :: i
        type(derivative_data_type_df5), dimension(5) :: uL_ddt, uR_ddt

        jac_L_R : do i = 1,2
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
                call roe_ddt(uL_ddt,uR_ddt,njk, dummy5,dfndu,wsn)
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

end module module_flux_jac_interface