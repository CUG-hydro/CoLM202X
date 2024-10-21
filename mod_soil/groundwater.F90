! 定义 find_jwt 函数
integer function find_jwt(zwt, zi_soisno, nl_soil)
  implicit none
  real(8), intent(in) :: zwt
  real(8), dimension(:), intent(in) :: zi_soisno(0:nl_soil)
  integer, intent(in) :: nl_soil
  integer :: j

  find_jwt = nl_soil
  ! allow jwt to equal zero when zwt is in top layer
  do j = 1, nl_soil
    if (zwt <= zi_soisno(j)) then
      find_jwt = j - 1
      exit
    endif
  end do
end function find_jwt


SUBROUTINE groundwater (nl_soil,deltim,pondmx,&
   eff_porosity,icefrac,&
   dz_soisno,zi_soisno,wice_soisno,wliq_soisno,&
   porsl,psi0,bsw,zwt,wa,&
   qcharge,rsubst,errw_rsub)

   USE MOD_Precision
   USE MOD_Const_Physical, only : tfrz

   ! ARGUMENTS:
   IMPLICIT NONE

   integer, intent(in) :: nl_soil       !
   real(r8), intent(in) :: deltim       ! land model time step (sec)
   real(r8), intent(in) :: pondmx       !

   real(r8), intent(in) :: eff_porosity(1:nl_soil)   ! effective porosity = porosity - vol_ice
   real(r8), intent(in) :: icefrac(1:nl_soil)        ! ice fraction (-)

   real(r8), intent(in) :: dz_soisno  (1:nl_soil)    ! layer depth (m)
   real(r8), intent(in) :: zi_soisno  (0:nl_soil)    ! interface level below a "z" level (m)
   real(r8), intent(inout) :: wice_soisno(1:nl_soil) ! ice lens (kg/m2)
   real(r8), intent(inout) :: wliq_soisno(1:nl_soil) ! liquid water (kg/m2)

   real(r8), intent(in) :: porsl(1:nl_soil)          ! volumetric soil water at saturation (porosity)
   real(r8), intent(in) :: psi0(1:nl_soil)           ! minimum soil suction (mm) [-]
   real(r8), intent(in) :: bsw(1:nl_soil)            ! Clapp and Hornberger "b"

   real(r8), intent(inout) :: zwt       ! the depth from ground (soil) surface to water table [m]
   real(r8), intent(inout) :: wa        ! water in the unconfined aquifer (mm)
   real(r8), intent(in)    :: qcharge   ! aquifer recharge rate (positive to aquifer) (mm/s)
   real(r8), intent(inout) :: rsubst    ! subsurface runoff (positive = out of soil column) (mm H2O /s)
   real(r8), intent(out)   :: errw_rsub ! the possible subsurface runoff dificit after PHS is included

   ! LOCAL ARGUMENTS

   integer  :: j                ! indices
   integer  :: jwt              ! index of the soil layer right above the water table (-)
   real(r8) :: xs               ! water needed to bring soil moisture to watmin (mm)
   real(r8) :: dzmm(1:nl_soil)  ! layer thickness (mm)
   real(r8) :: xsi              ! excess soil water above saturation at layer i (mm)
   real(r8) :: xsia             ! available pore space at layer i (mm)
   real(r8) :: xs1              ! excess soil water above saturation at layer 1 (mm)
   real(r8) :: ws               ! summation of pore space of layers below water table (mm)
   real(r8) :: s_node           ! soil wetness (-)
   real(r8) :: available_wliq_soisno     ! available soil liquid water in a layer
   real(r8) :: qcharge_tot      !
   real(r8) :: qcharge_layer    !
   real(r8) :: drainage         !
   real(r8) :: drainage_tot     !
   real(r8) :: drainage_layer   !
   real(r8) :: s_y              !
   real(r8) :: sy               ! specific yield [-]

   real(r8) :: wt
   real(r8) :: wtsub
   real(r8) :: dzsum
   real(r8) :: icefracsum
   real(r8) :: fracice_rsub
   real(r8) :: imped

   real(r8), parameter :: watmin = 0.01  ! Limit irreduciable wrapping liquid water
   ! a tunable constant
   real(r8), parameter :: rsbmx  = 5.0   ! baseflow coefficient [mm/s]
   real(r8), parameter :: timean = 10.5  ! global mean topographic index


   ! -------------------------------------------------------------------------
   !   ! Convert layer thicknesses from m to mm
   DO j = 1,nl_soil
      dzmm(j) = dz_soisno(j)*1000.
   ENDDO

   !     ! The layer index of the first unsaturated layer,
   !     ! i.e., the layer right above the water table
   jwt = nl_soil
   ! allow jwt to equal zero when zwt is in top layer
   DO j = 1, nl_soil
      IF(zwt <= zi_soisno(j)) THEN
         jwt = j-1
         EXIT
      ENDIF
   ENDDO

   !============================== QCHARGE =========================================
   ! Water table changes due to qcharge
   ! use analytical expression for aquifer specific yield
   sy = porsl(nl_soil)*(1.-(1.-1.e3*zwt/psi0(nl_soil))**(-1./bsw(nl_soil)))
   sy = max(sy,0.02)

   wa = wa + qcharge*deltim
   !
   !---------------------------------------
   ! water table is below the soil column
   IF(jwt == nl_soil) THEN
      zwt = max(0.,zwt - (qcharge*deltim)/1000./sy)
   ELSE
      ! water table within soil layers 1-9
      ! try to raise water table to account for qcharge
      qcharge_tot = qcharge * deltim

      IF(qcharge_tot > 0.) THEN ! rising water table
         DO j = jwt+1, 1,-1
            ! use analytical expression for specific yield

            s_y = porsl(j) * (1.-(1.-1.e3*zwt/psi0(j))**(-1./bsw(j)))
            s_y=max(s_y,0.02)

            qcharge_layer = min(qcharge_tot,(s_y*(zwt-zi_soisno(j-1))*1.e3))
            qcharge_layer = max(qcharge_layer,0.)

            zwt = max(0.,zwt - qcharge_layer/s_y/1000.)

            qcharge_tot = qcharge_tot - qcharge_layer
            IF (qcharge_tot <= 0.) EXIT
         ENDDO
      ELSE ! deepening water table (negative qcharge)
         DO j = jwt+1, nl_soil
            ! use analytical expression for specific yield
            s_y = porsl(j) * (1.-(1.-1.e3*zwt/psi0(j))**(-1./bsw(j)))
            s_y=max(s_y,0.02)
            qcharge_layer = max(qcharge_tot,-(s_y*(zi_soisno(j) - zwt)*1.e3))
            qcharge_layer = min(qcharge_layer,0.)
            qcharge_tot = qcharge_tot - qcharge_layer ! 水位下降则出水

            IF (qcharge_tot >= 0.) THEN
               zwt = max(0.,zwt - qcharge_layer/s_y/1000.)
               EXIT
            ELSE
               zwt = zi_soisno(j)
            ENDIF
         ENDDO
         IF (qcharge_tot > 0.) zwt = max(0.,zwt - qcharge_tot/1000./sy)
      ENDIF
   ENDIF

   !-- Topographic runoff  ----------------------------------------------------------
   IF (DEF_Runoff_SCHEME == 0) THEN
      CALL SubsurfaceRunoff_SIMTOP (nl_soil, icefrac, dz_soisno, zi_soisno, zwt, rsubst)
   ENDIF

   drainage = rsubst
   ! dzsum = 0.
   ! icefracsum = 0.
   ! DO j = max(jwt,1), nl_soil
   !    dzsum = dzsum + dzmm(j)
   !    icefracsum = icefracsum + icefrac(j) * dzmm(j)
   ! ENDDO
   ! ! add ice impedance factor to baseflow
   ! fracice_rsub = max(0.,exp(-3.*(1.-(icefracsum/dzsum)))-exp(-3.))/(1.0-exp(-3.))
   ! imped = max(0.,1.-fracice_rsub)
   ! drainage = imped * 5.5e-3 * exp(-2.5*zwt)  ! drainage (positive = out of soil column)

   !-- Water table is below the soil column  ----------------------------------------
   IF(jwt == nl_soil) THEN
      wa = wa - drainage * deltim
      zwt = max(0.,zwt + (drainage * deltim)/1000./sy)
      wliq_soisno(nl_soil) = wliq_soisno(nl_soil) + max(0.,(wa-5000.))
      wa = min(wa, 5000.) ! limit water in aquifer to 5000 mm
   ELSE
      !-- Water table within soil layers 1-9  ------------------------------------------
      !============================== RSUB_TOP =========================================
      !-- Now remove water via drainage
      drainage_tot = - drainage * deltim
      DO j = jwt+1, nl_soil
         ! use analytical expression for specific yield
         s_y = porsl(j) * ( 1. - (1.-1.e3*zwt/psi0(j))**(-1./bsw(j)))
         s_y = max(s_y,0.02)

         drainage_layer = max(drainage_tot, -(s_y*(zi_soisno(j)-zwt)*1.e3))
         drainage_layer = min(drainage_layer,0.)
         wliq_soisno(j) = wliq_soisno(j) + drainage_layer

         drainage_tot = drainage_tot - drainage_layer

         IF(drainage_tot >= 0.)THEN
            zwt = max(0.,zwt - drainage_layer/s_y/1000.)
            EXIT
         ELSE
            zwt = zi_soisno(j)
         ENDIF
      ENDDO

      !-- Remove residual drainage  ------------------------------------------------
      zwt = max(0.,zwt - drainage_tot/1000./sy)
      wa = wa + drainage_tot

      !-- Recompute jwt  ---------------------------------------------------------------
      ! allow jwt to equal zero when zwt is in top layer
      jwt = nl_soil
      DO j = 1, nl_soil
         IF(zwt <= zi_soisno(j)) THEN
            jwt = j-1
            EXIT
         ENDIF
      ENDDO

   ENDIF   ! end of jwt IF construct

   zwt = max(0.0,zwt)
   zwt = min(80.,zwt)

   rsubst = drainage

   ! Correction [1]
   ! NON-physically based corection on wliq_soisno
   ! excessive water above saturation added to the above unsaturated layer like a bucket
   ! IF column over saturated, excess water goes to runoff

   DO j = nl_soil,2,-1
      xsi = max(wliq_soisno(j)-eff_porosity(j)*dzmm(j),0.)
      wliq_soisno(j) = min(eff_porosity(j)*dzmm(j), wliq_soisno(j))
      wliq_soisno(j-1) = wliq_soisno(j-1) + xsi
   ENDDO

   ! 12/2022, note by yuan: a potential bug below which needs check,
   ! if wice_soisno(1) > pondmx + porsl*dzmm, so xs1>0, in that case,
   ! wliq_soisno(1) will be nagtive, and xs1 is positive.
   xs1 = wliq_soisno(1) - (pondmx+porsl(1)*dzmm(1)-wice_soisno(1))
   IF(xs1 > 0.)THEN
      wliq_soisno(1) = pondmx+porsl(1)*dzmm(1)-wice_soisno(1)
   ELSE
      xs1 = 0.
   ENDIF

   rsubst = rsubst + xs1 / deltim
   ! Correction [2]
   ! NON-physically based corection on wliq_soisno
   ! Limit wliq_soisno to be greater than or equal to watmin.
   ! Get water needed to bring wliq_soisno equal watmin from lower layer.
   ! If insufficient water in soil layers, get from aquifer water
   xs = 0.
   DO j = 1, nl_soil
      IF (wliq_soisno(j) < 0.) THEN
         xs = xs + wliq_soisno(j)
         wliq_soisno(j) = 0.
      ENDIF
   ENDDO

   ! Sub-surface runoff and drainage
   errw_rsub = min(0., rsubst + xs/deltim)
   rsubst = max(0., rsubst + xs/deltim)
END SUBROUTINE groundwater
