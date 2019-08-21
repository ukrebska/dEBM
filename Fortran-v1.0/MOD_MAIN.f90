MODULE MOD_MAIN
! ************************************************************************
! * This MODULE calculates surface mass balance                          *
! ************************************************************************
  USE MOD_PRE
  USE MOD_DATA

  IMPLICIT NONE

contains

SUBROUTINE CALC_ANGLES(winkelns, winkelds, winkelws, winkeli)
  real, intent(out) :: winkelns, winkelds, winkelws, winkeli
  ! angles
  winkelns = asin( -c2cs /( 1 - Ans )/( 1330 * taucs ))* 180/pi
  winkelds = asin( -c2cs /( 1 - Ads )/( 1330 * taucs ))* 180/pi
  winkelws = asin( -c2cs /( 1 - Aws )/( 1330 * taucs ))* 180/pi
  winkeli  = asin( -c2cs /( 1 - Ai  )/( 1330 * taucs ))* 180/pi
END SUBROUTINE CALC_ANGLES


SUBROUTINE dEBM_decl(days_in, obl_in, DECL_out)
  ! ************************************************************************
  ! * so far this is a very simple calculation                             *
  ! * which assumes eccentricity is 1...                                   *
  ! * and spring equinox  is on March 20th                                 *
  ! * Values get more inaccurate towards autumn, lap years are neglected   *
  ! * In the present day orbital configuration, the Earth orbits slower    *
  ! * in boreal summer and decl=23.44*sin(pi/180*186/180*(days(:)-79)) is  *
  ! * more accurate for the northern hemisphere but not a good choice      *
  ! * for the southern hemisphere...                                       *
  ! * uncomment the following to see the difference                        *
  ! ************************************************************************
  ! * INPUT: days, obliqu                                                  *
  ! *   days           :  days in each month                               *
  ! *   obliqu         :  obliquity                                        *
  ! ************************************************************************
  ! * OUTPUT: DECL                                                         *
  ! *   DECL          : declinity                                          *
  ! ************************************************************************

  real, intent(in) :: days_in
  real, intent(in) :: obl_in
  real, intent(out) :: DECL_out

  DECL_out = obl_in * sin(pi/180.*360./365.*(days_in-79.))   ! 79 is the number of days since New Year for the spring equinox
  ! write(*,*) "obl_in", obl_in
  ! write(*,*) "days_in", days_in
  ! write(*,*) "sin(pi/180.*360./365.*(days_in-79.))", sin(pi/180.*360./365.*(days_in-79.))

END SUBROUTINE dEBM_decl


SUBROUTINE dEBM_sunny_hours_c(mth_in, elev_in, lat_in, obl_in, HOURS_OUT, Q_OUT, FLUXFAC_OUT)
  ! ************************************************************************
  ! calculates the time, the sun is above a certain elevation angle        *
  ! ************************************************************************
  ! * INPUT: lat, elev                                                     *
  ! *   lat           : latitude                                           *
  ! *   elev          : elevation                                          *
  ! ************************************************************************
  ! * OUTPUT: HOURS, Q, FLUXFAC                                            *
  ! *   HOURS         : duration of melt period in hours                   *
  ! *   Q             : daily insolation                                   *
  ! *   FLUXFAC       : daily insolation SW                                *
  ! ************************************************************************

  real, intent(in) :: elev_in, obl_in
  real, intent(in), dimension(:,:) :: lat_in
  integer, intent(in) :: mth_in

  integer :: i, j, k

  real, allocatable, dimension(:,:) :: elev1, decl
  real, allocatable, dimension(:,:) :: sinphisind,cosphicosd
  real, allocatable, dimension(:,:) :: ha_0,ha_elev
  complex, allocatable, dimension(:,:) :: ha_elev0
  real, allocatable, dimension(:,:) :: sol_flux_fact_0,sol_flux_fact_e

  real, intent(out), dimension(:,:) :: HOURS_OUT, Q_OUT, FLUXFAC_OUT

  allocate (decl(xlen,ylen))
  allocate (elev1(xlen,ylen))
  allocate (sinphisind(xlen,ylen))
  allocate (cosphicosd(xlen,ylen))
  allocate (ha_0(xlen,ylen))
  allocate (ha_elev(xlen,ylen))
  allocate (ha_elev0(xlen,ylen))
  allocate (sol_flux_fact_0(xlen,ylen))
  allocate (sol_flux_fact_e(xlen,ylen))

  ! decl
  CALL dEBM_decl(days(mth_in), obl_in, decl(1,1))

  do i = 1, xlen
      do j = 1, ylen
          decl(i, j) = decl(1,1)
      end do
  end do

! write(*,*) "decl", decl(1,1)

  ! elev
  elev1(:,:)       = elev_in*pi/180.
! write(*,*) "elev", elev1(448, 155)

! ! elev new
  ! elev1 = (decl + (90. - lat_in))*pi/180.
!   ! write(*,*) "elev1", elev1(448, 155)
!   ! write(*,*) "elev_in", elev_in

  ! !products of sin(lat)*sin(decl) and cos(lat)*cos(decl)  is needed more than once
  sinphisind       = sin(pi/180*lat_in)*sin(pi/180*decl)
  cosphicosd       = cos(pi/180*lat_in)*cos(pi/180*decl)
  ! write(*,*) "sinphisind", sinphisind(448, 155)
  ! write(*,*) "cosphicosd", cosphicosd(448, 155)

  ! hourangle of sun rise and factor which relates solar flux density and
  ! daily insolation SW
   ! ha_0             = real(acos(-sinphisind/cosphicosd))
   ! sol_flux_fact_0  = (sinphisind*ha_0+cosphicosd*sin(ha_0))/pi
   ! write(*,*) "ha_0", ha_0(448, 155)
   ! write(*,*) "sol_flux_fact_0", sol_flux_fact_0(448, 155)

   ! ha_0 = 0.
   ! sol_flux_fact_0  = 0.
   ! where (abs(sinphisind/cosphicosd)<1.)
   ha_0 = -sinphisind/cosphicosd
     where (ha_0 > 1.) ha_0 = 1.
     where (ha_0 < -1.) ha_0 = -1.
     ha_0 = acos(ha_0)
  sol_flux_fact_0  = (sinphisind*ha_0+cosphicosd*sin(ha_0))/pi
   ! end where
   ! write(*,*) "-sinphisind/cosphicosd", -sinphisind(448, 155)/cosphicosd(448, 155)
   ! write(*,*) "ha_0", ha_0(448, 155)
   ! write(*,*) "sol_flux_fact_0", sol_flux_fact_0(448, 155)


  ! hourangle of sun rising above elev
  ha_elev = (sin(elev1)-sinphisind)/cosphicosd
    where (ha_elev > 1.) ha_elev = 1.
    where (ha_elev < -1.) ha_elev = -1.
  ha_elev = acos(ha_elev)
    ! write(*,*) "(sin(elev1))", (sin(elev1))
    ! write(*,*) "(sin(elev1)-sinphisind)", (sin(elev1)-sinphisind(448, 155))
    ! write(*,*) "(sin(elev1)-sinphisind)/cosphicosd", (sin(elev1)-sinphisind(448, 155))/cosphicosd(448, 155)
    ! write(*,*) "acos((sin(elev1)-sinphisind)/cosphicosd)", acos((sin(elev1)-sinphisind(448, 155))/cosphicosd(448, 155))
    ! ha_elev(:,:)           = acos(0.8)
    ! write(*,*) "ha_elev0", ha_elev0(448, 155)
    ! write(*,*) "ha_elev", ha_elev(448, 155)
    ! write(*,*) "ha_elev", maxval(ha_elev)
    ! write(*,*) "ha_elev", minval(ha_elev)


  ! duration of melt period in hours
  HOURS_OUT         = ha_elev*24/pi
  ! write(*,*) "HOURS_OUT", HOURS_OUT(448, 155)
  ! write(*,*) "HOURS_OUT", minval(HOURS_OUT)

  ! proportion of daily insolation which is effective during melt period
  sol_flux_fact_e   = (sinphisind*ha_elev+cosphicosd*sin(ha_elev))/pi
  ! Q_OUT             = sol_flux_fact_e/sol_flux_fact_0

  where (sol_flux_fact_0==0.)
    Q_OUT = 0.
  elsewhere
    Q_OUT = sol_flux_fact_e/sol_flux_fact_0
  end where

  ! write(*,*) "sol_flux_fact_e", sol_flux_fact_e(448, 155)
  ! write(*,*) "sol_flux_fact_0", sol_flux_fact_0(448, 155)
  ! write(*,*) "Q_OUT", Q_OUT(448, 155)

  FLUXFAC_OUT       = sol_flux_fact_0
  ! write(*,*) "FLUXFAC_OUT", FLUXFAC_OUT(448, 155)

END SUBROUTINE dEBM_sunny_hours_c



SUBROUTINE PDD4(temp_in, stddev_in, PDD_OUT)
  ! ************************************************************************
  ! * calculates the time, the sun is above a certain elevation angle      *
  ! * PDD is approximated as described in                                  *
  ! * Calov R and Greve R (2005)                                           *
  ! * Correspondence. A semi-analytical solution for the positive          *
  ! * degree-day model with stochastic temperature variations.             *
  ! * J. Glaciol., 51 (172), 173?175                                       *
  ! * (doi:10.3189/172756505781829601)                                     *
  ! ************************************************************************
  ! * INPUT: Temp, stddev                                                  *
  ! *   Temp          : surface temperature                                *
  ! *   stddev        : elevation                                          *
  ! ************************************************************************
  ! * OUTPUT: PDD                                                          *
  ! *   PDD           : duration of melt period in hours                   *
  ! ************************************************************************

  real, intent(in) :: stddev_in
  real, intent(in), dimension(:,:) :: temp_in
  real, intent(out), dimension(:,:) :: PDD_OUT

  PDD_OUT = stddev_in/sqrt(2*pi)*exp(-(temp_in**2)/(2*(stddev_in**2)))+temp_in/2*erfc(-temp_in/(sqrt(2.)*stddev_in))

END SUBROUTINE PDD4



SUBROUTINE dEBMmodel_cloud(A_in, Aoc_in, swd_in, Temp_in, PDD_in, rain_in, hours_in, q_in, c1cs_in, c2cs_in, c1oc_in, c2oc_in, MC_in, fluxfac_in, taucs_in, tauoc_in, MELT_OUT, REFR_OUT)
! SUBROUTINE dEBMmodel_cloud(A_in, Aoc_in, swd_in, Temp_in, PDD_in, rain_in, hours_in, q_in, fluxfac_in, taucs_in, tauoc_in, MELT_OUT, REFR_OUT, LWD_EST_OUT,TAU_OUT)
  ! ************************************************************************
  ! * calculates surface melt and refreeze                                 *
  ! ************************************************************************
  ! * INPUT: A,swd,T,PDD,rain,hours,q,MC,fluxfac                           *
  ! *   A             : Albedo                                             *
  ! *   Aoc             : Albedo                                             *
! *   swd           : daily short wave insolation at the surfacet        *
  ! *   PDD           : an approx. of melt-period temperature              *
  ! *   rain          : liquid precipitation                               *
  ! *   hours         : duration of melt period in hours                   *
  ! *   q             : daily insolation                                   *
  ! *   MC            : (logic) the background melt condition (e.g. T>-6.5)*
  ! *   fluxfac       : daily insolation SW                                *
! **************************************************************************
  ! * OUTPUT: MELT, REFR                                                   *
  ! *   MELT          : the melt rate (mm/day)                             *
  ! *   REFR          : refreeze                                           *
  ! *   LWD_EST_OUT          : refreeze                                           *
  ! *   TAU_OUT          : refreeze                                           *
  ! ************************************************************************

  real, intent(in) :: A_in, Aoc_in, taucs_in, tauoc_in
  real, intent(in) :: c1cs_in, c2cs_in, c1oc_in, c2oc_in
  real, intent(in), dimension(:,:)    :: swd_in, Temp_in, PDD_in, rain_in, hours_in, q_in, fluxfac_in
  logical, intent(in), dimension(:,:) :: MC_in
  real, intent(out), dimension(:,:)   :: MELT_OUT, REFR_OUT
 ! real, intent(out), dimension(:,:)   ::  LWD_EST_OUT,TAU_OUT

  real, parameter :: Lf  = 3.34e5         ! Lf is latent heat of fusion
  real, parameter :: S0  = 1330.                    ! solar constant
  real, allocatable, dimension(:,:)   :: Scs, Soc, Scsr, Socr, gwi
  real, allocatable, dimension(:,:)   :: meltcs, meltoc, refroc, refrcs, meltFDcs, meltFDoc

  allocate (Scs(xlen,ylen))
  allocate (Soc(xlen,ylen))
  allocate (Scsr(xlen,ylen))
  allocate (Socr(xlen,ylen))
  allocate (gwi(xlen,ylen))
  allocate (meltcs(xlen,ylen))
  allocate (meltoc(xlen,ylen))
  allocate (refrcs(xlen,ylen))
  allocate (refroc(xlen,ylen))
  allocate (meltFDcs(xlen,ylen))
  allocate (meltFDoc(xlen,ylen))

  Scs = fluxfac_in*S0*taucs_in        ! expected sunshine for clear sky conditions
  Soc = fluxfac_in*S0*tauoc_in        ! expected sunshine for completely overcast conditions
  ! write(*,*) "fluxfac_in", fluxfac_in(448, 155)
  ! write(*,*) "Scs", Scs(448, 155)
  ! write(*,*) "Soc", Soc(448, 155)
  ! ?????
  ! why use Scsr, the maximum of real and expected sunshine, rather than real
  Scsr = max(Scs,swd_in)
  Socr = min(Soc,swd_in);
  ! ! write(*,*) "swd_in", swd_in(448, 155)
  ! write(*,*) "Scsr", Scsr(448, 155)
  ! write(*,*) "Socr", Socr(448, 155)

  ! good weather index
  ! percentage of days in a month, which has clear sky conditions
  ! (there is just overcast or clear conditions, nothing inbetween)
  gwi = Scs*0. + 1.
  ! TAU_OUT = Scs*asin(1.2)
  where (swd_in>50.)
    ! TAU_OUT = swd_in./(fluxfac_in*S0);
    gwi = max(0., min(1., (swd_in - Soc)/(Scs - Soc)))
  end where
  ! write(*,*) "swd_in", swd_in(448, 155)
  ! ! write(*,*) "swd_in - Soc", (swd_in(448, 155) - Soc(448, 155))
  ! ! write(*,*) "Scs - Soc", (Scs(448, 155) - Soc(448, 155))
  ! write(*,*) "gwi", gwi(448, 155)
!
  ! 0.05 is penetration
  meltcs = 0.
  where (MC_in) meltcs   = (((1-A_in)*Scsr*24.*q_in)+(c2cs+c1cs*PDD_in)*hours_in)/Lf*60.*60.
  ! write(*,*) "A_in", A_in
  ! write(*,*) "Scsr", Scsr(448, 155)
  ! write(*,*) "q_in", q_in(448, 155)
  ! write(*,*) "PDD_in", PDD_in(448, 155)
  ! write(*,*) "hours_in", hours_in(448, 155)
  ! ! write(*,*) "MC_in", MC_in(448, 155)
  ! write(*,*) "meltcs", meltcs(448, 155)

  ! meltFDcs
  meltFDcs = 0.
  where (MC_in) meltFDcs = ((1-A_in)*Scsr+(c2cs+c1cs*Temp_in))/Lf*24.*60.*60.
  ! meltFDoc
  meltFDoc = 0.
  where (MC_in) meltFDoc = ((1-Aoc_in)*Socr+(c2oc+(c1oc)*Temp_in))/Lf*24.*60.*60.        ! we don't distinguish night and day here...

  refroc = max(-meltFDoc, 0.)
  meltoc = max( meltFDoc, 0.)
  meltcs = max(0., max(meltFDcs, meltcs))
  refrcs = max(0., meltcs - meltFDcs)

  ! MELT_OUT
  MELT_OUT = gwi*meltcs + (1.-gwi)*meltoc;
  ! REFR_OUT
  REFR_OUT = min(MELT_OUT + rain_in, (gwi*refrcs + (1.-gwi)*refroc))
  ! ! LWD_EST_OUT
  ! LWD_EST_OUT = bolz * (gwi*epsacs + (1-gwi)*epsaoc) * ((T+T0)**4.)

  ! ! MELT_cs
  ! where (isnan(meltcs))
  !   MELT_cs = 0.
  ! elsewhere
  !   MELT_cs = gwi*max(0.,meltcs)
  ! end where
  ! ! MELT_oc
  ! where (isnan(gwi*(meltFDcs-meltcs)+(1-gwi)*meltFDoc))
  !   MELT_oc = 0.
  ! elsewhere
  !   MELT_oc = max(0.,gwi*(meltFDcs-meltcs)+(1-gwi)*meltFDoc)
  ! end where
  ! !
  ! MELT_OUT = MELT_cs + MELT_oc
  ! write(*,*) "MELT_cs", MELT_cs(448, 155)
  ! write(*,*) "MELT_oc", MELT_oc(448, 155)

  ! write(*,*) "MC_in", MC_in(448, 155)
  ! write(*,*) "meltcs", meltcs(448, 155)
  ! write(*,*) "meltFDcs", meltFDcs(448, 155)
  ! write(*,*) "meltFDoc", meltFDoc(448, 155)
  ! write(*,*) "MELT_OUT", MELT_OUT(448, 155)

  ! REFR_OUT = max(0.,min(MELT_OUT + rain_in, gwi*(meltcs-meltFDcs)))
  ! where (isnan(gwi*(meltcs-meltFDcs)))
  !   REFR_OUT = MELT_OUT + rain_in
  ! elsewhere (isnan(MELT_OUT + rain_in))
  !   REFR_OUT = gwi*(meltcs-meltFDcs)
  ! elsewhere
  !   REFR_OUT = min(MELT_OUT + rain_in, gwi*(meltcs-meltFDcs))
  ! end where
  ! !
  ! where (isnan(REFR_OUT))
  !   REFR_OUT = 0.
  ! elsewhere
  !   REFR_OUT = max(0.,REFR_OUT)
  ! end where

  ! write(*,*) "REFR_OUT", REFR_OUT(448, 155)

END SUBROUTINE dEBMmodel_cloud




SUBROUTINE dEBM_core(tempm_in, radm_in, ppm_in, tmpSNH_in, lastyear_in, lat_in, mask_in, obl_in, SNH_OUT, SMB_OUT, MELT_OUT, ACC_OUT, REFR_OUT, A_OUT)
  ! ************************************************************************
  ! * This SUBTROUTINE calculates the surface mass blance                  *
  ! ************************************************************************
  ! * INPUT: tempm, radm, ppm, tmpSNH, lastyear, lat, msk                  *
  ! *  tempm,radm,ppm are 3-dimension                                      *
  ! *  tmpSNH, lastyear,lat,msk are 2 dimension                            *
  ! *   tempm         : temperature degC                                       *
  ! *   radm          : daily short wave insolation at the surface         *
  ! *   ppm           : precipitation mm/day                                      *
  ! *   tmpSNH        : snow height of the last timestep                   *
  ! *   lastyear      : snow height at the end of last summer              *
  ! *   lat           : latitude                                           *
  ! *   mask          : mapping                                            *
  ! *   obl          : obliquity                                            *
  ! *   lhf          :                                             *
  ! ************************************************************************
  ! * OUTPUT: SNH, SMB, MELT, ACC, REFR, A                                 *
  ! *   SNH           : snow height                                        *
  ! *   SMB           : surface mass balance                               *
  ! *   MELT          : surface melt                                       *
  ! *   ACC           : accumulation                                       *
  ! *   REFR          : refreeze                                           *
  ! *   A             : albedo                                             *
  ! ************************************************************************

  real, intent(in), dimension(:,:,:) :: tempm_in, radm_in, ppm_in, lat_in
  real, intent(in), dimension(:,:)   :: tmpSNH_in, lastyear_in
  logical, intent(in), dimension(:,:,:) :: mask_in
  real, intent(in) :: obl_in
  real, intent(out), dimension(:,:,:) :: SNH_OUT, SMB_OUT, MELT_OUT, ACC_OUT, REFR_OUT, A_OUT

  integer :: k
  real, allocatable, dimension(:,:) :: hoursns, qns, fluxfacns,&
                                          &hoursds, qds, fluxfacds,&
                                          &hoursws, qws, fluxfacws,&
                                          &hoursi, qi, fluxfaci,&
                                          &PDD0
  real, allocatable, dimension(:,:) :: solid0, solid, rain, snow
  logical, allocatable, dimension(:,:) :: tempmask
  real, allocatable, dimension(:,:) :: MELTns, PREFRns
  real, allocatable, dimension(:,:) :: MELTds, PREFRds
  real, allocatable, dimension(:,:) :: MELTws, PREFRws
  real, allocatable, dimension(:,:) :: MELTi, PREFRi
  logical, allocatable, dimension(:,:) :: old_wet, new_snow, dry_snow, wet_snow
  real, allocatable, dimension(:,:) :: melt_snow, ice_melt, ice_melt_fract
  real, allocatable, dimension(:,:) :: Asnow
  real, allocatable, dimension(:,:)   :: tmpSNH, lastyear

  !
  allocate (hoursns(xlen,ylen))
  allocate (qns(xlen,ylen))
  allocate (fluxfacns(xlen,ylen))
  !
  allocate (hoursds(xlen,ylen))
  allocate (qds(xlen,ylen))
  allocate (fluxfacds(xlen,ylen))
  !
  allocate (hoursws(xlen,ylen))
  allocate (qws(xlen,ylen))
  allocate (fluxfacws(xlen,ylen))
  !
  allocate (hoursi(xlen,ylen))
  allocate (qi(xlen,ylen))
  allocate (fluxfaci(xlen,ylen))
  !
  allocate (PDD0(xlen,ylen))
  !
  allocate (solid0(xlen,ylen))
  allocate (solid(xlen,ylen))
  allocate (rain(xlen,ylen))
  allocate (snow(xlen,ylen))
  !
  allocate (tempmask(xlen,ylen))
  !
  allocate (MELTns(xlen,ylen))
  allocate (PREFRns(xlen,ylen))
  allocate (MELTds(xlen,ylen))
  allocate (PREFRds(xlen,ylen))
  allocate (MELTws(xlen,ylen))
  allocate (PREFRws(xlen,ylen))
  allocate (MELTi(xlen,ylen))
  allocate (PREFRi(xlen,ylen))
  !
  allocate (old_wet(xlen,ylen))
  allocate (new_snow(xlen,ylen))
  allocate (dry_snow(xlen,ylen))
  allocate (wet_snow(xlen,ylen))
  allocate (melt_snow(xlen,ylen))
  !
  allocate (ice_melt(xlen,ylen))
  allocate (ice_melt_fract(xlen,ylen))
  allocate (Asnow(xlen,ylen))
  allocate (tmpSNH(xlen,ylen))
  allocate (lastyear(xlen,ylen))

  ! write(*,*) "calculating parameters c1 anc c2..."
  CALL CALC_PARAMETER(c1cs, c2cs, c1oc, c2oc)
  ! write(*,*) "c1cs, c2cs, c1oc, c2oc = ",c1cs, c2cs, c1oc, c2oc

  ! write(*,*) "calculating angles...."
  CALL CALC_ANGLES(winkelns, winkelds, winkelws, winkeli)
  ! write(*,*) "winkelns, winkelds, winkelws, winkeli",winkelns, winkelds, winkelws, winkeli

  !
  SMB_OUT    = tempm_in*0.
  SNH_OUT    = tempm_in*0.
  MELT_OUT   = tempm_in*0.
  ACC_OUT    = tempm_in*0.
  REFR_OUT   = tempm_in*0.
  A_OUT      = tempm_in*0.
  wet_snow(:,:) = .FALSE.

  tmpSNH = tmpSNH_in
  lastyear = lastyear_in

  do k = 1, 12

    ! write(*,*) "calculates month",k

    ! write(*,*) "calculates the time that the sun is above a certain elevation angle...."
    CALL dEBM_sunny_hours_c(k, winkelns, lat_in(:,:,k), obliquity, hoursns, qns, fluxfacns)
    CALL dEBM_sunny_hours_c(k, winkelds, lat_in(:,:,k), obliquity, hoursds, qds, fluxfacds)
    CALL dEBM_sunny_hours_c(k, winkelws, lat_in(:,:,k), obliquity, hoursws, qws, fluxfacws)
    CALL dEBM_sunny_hours_c(k, winkeli, lat_in(:,:,k), obliquity, hoursi, qi, fluxfaci)

    ! write(*,*) k, " hoursns", hoursns(448, 155,k)
    ! write(*,*) k, " qns ", qns(448, 155,k)
    ! write(*,*) k, " fluxfacns ", fluxfacns(448, 155,k)

    ! to determine the fraction of solid and liquid water in precipitation
    solid0 = 0.
    where (tempm_in(:,:,k)>-slim) solid0 = sin(tempm_in(:,:,k)/2./slim*pi)+1.
    solid = 0.
    where (tempm_in(:,:,k)<slim) solid = 1.-.5*solid0
    ! solid(:,:,k) = max(0.,solid(:,:,k))
    ! write(*,*) "tempm_in",tempm_in(448, 155,k)
    ! write(*,*) "slim",slim
    ! write(*,*) "solid0",solid0(448, 155,k)
    ! write(*,*) "solid",solid(448, 155,k)
    ! write(*,*) "ppm_in",ppm_in(448, 155,k)

    snow = 0.
    where (mask_in(:,:,k)) snow  = solid * ppm_in(:,:,k)            ! snow&rain  fractionation acc. to Pipes & Quick 1977
    ! write(*,*) "snow",snow(448, 155,k)
    rain = 0.
    where (mask_in(:,:,k)) rain  = (1 - solid) * ppm_in(:,:,k)
    ! write(*,*) "rain",rain(448, 155,k)

    tempmask = ((tempm_in(:,:,k) > -6.5) .AND. (mask_in(:,:,k)))
    ! write(*,*) "temp",tempm_in(448, 155,k)
    ! write(*,*) "ppm",ppm_in(448, 155,k)
    ! write(*,*) "radm_in",radm_in(448, 155,k)
    ! write(*,*) "snow",snow(448, 155,k)
    ! write(*,*) "rain",rain(448, 155,k)
    ! write(*,*) "mask_in",mask_in(448, 155,k)
    ! write(*,*) "tempmask",tempmask(448, 155,k)

    ! PDD=PDD4(temp,3.5)*msk
    ! write(*,*) "calculates PDD...."
    CALL PDD4(tempm_in(:,:,k), stddev, PDD0)
    ! write(*,*) "PDD0",PDD0(448, 155,k)

    ! write(*,*) "calculates dEBMmodel_cloud...."
    ! write(*,*) "Ans"
    CALL dEBMmodel_cloud(Ans, Ans+.05, radm_in(:,:,k), tempm_in(:,:,k), PDD0, rain, hoursns, qns, c1cs, c2cs, c1oc, c2oc, tempmask, fluxfacns, taucs, tauoc,  MELTns, PREFRns)
    ! write(*,*) "Ans",Ans
    ! write(*,*) "MELTns",MELTns(448, 155)
    ! write(*,*) "PREFRns",PREFRns(448, 155)
    CALL dEBMmodel_cloud(Ads, Ads+.05, radm_in(:,:,k), tempm_in(:,:,k), PDD0, rain, hoursds, qds, c1cs, c2cs, c1oc, c2oc, tempmask, fluxfacns, taucs, tauoc, MELTds, PREFRds)
    ! write(*,*) "Ads",Ads
    ! write(*,*) "MELTds",MELTds(448, 155)
    ! write(*,*) "PREFRds",PREFRds(448, 155)
    CALL dEBMmodel_cloud(Aws, Aws+.05, radm_in(:,:,k), tempm_in(:,:,k), PDD0, rain, hoursws, qws, c1cs, c2cs, c1oc, c2oc, tempmask, fluxfacns, taucs, tauoc, MELTws, PREFRws)
    ! write(*,*) "Aws",Aws
    ! write(*,*) "MELTws",MELTws(448, 155)
    ! write(*,*) "PREFRws",PREFRws(448, 155)
    CALL dEBMmodel_cloud(Ai , Ai+.05,  radm_in(:,:,k), tempm_in(:,:,k), PDD0, rain, hoursi,  qi, c1cs, c2cs, c1oc, c2oc, tempmask, fluxfacns, taucs, tauoc, MELTi, PREFRi)
    ! write(*,*) "Ai",Ai
    ! write(*,*) "MELTi",MELTi(448, 155)
    ! write(*,*) "PREFRi",PREFRi(448, 155)

    ! write(*,*) k, " hoursns", hoursns(448, 155,k)
    ! write(*,*) k, " hoursds", hoursds(448, 155,k)
    ! write(*,*) k, " hoursws", hoursws(448, 155,k)
    ! write(*,*) k, " hoursi", hoursi(448, 155,k)

    ! write(*,*) k, " fluxfacns", fluxfacns(448, 155,k)
    ! write(*,*) k, " fluxfacds", fluxfacds(448, 155,k)
    ! write(*,*) k, " fluxfacws", fluxfacws(448, 155,k)
    ! write(*,*) k, " fluxfaci", fluxfaci(448, 155,k)

    ! write(*,*) "calculates snow type...."
    old_wet  = (wet_snow .AND. (PREFRws < (MELTws + rain)))
    ! write(*,*) "wet_snow",wet_snow(448, 155,k)
    ! write(*,*) "PREFRws",PREFRws(448, 155,k)
    ! write(*,*) "MELTws",MELTws(448, 155,k)
    ! write(*,*) "rain",rain(448, 155,k)
    ! write(*,*) "old_wet",old_wet(448, 155,k)
    new_snow = (MELTns <= snow)
    dry_snow = ((.NOT.new_snow) .AND. (PREFRds >= (MELTds + rain)) .AND. (.NOT.old_wet))
    wet_snow = ((.NOT.new_snow) .AND. ((.NOT.dry_snow) .OR. (old_wet)))
    ! write(*,*) "MELTns",MELTns(448, 155,k)
    ! write(*,*) "snow",snow(448, 155,k)
    ! write(*,*) "rain",rain(448, 155,k)
    ! write(*,*) "old_wet",old_wet(448, 155)
    ! write(*,*) "new_snow",new_snow(448, 155)
    ! write(*,*) "dry_snow",dry_snow(448, 155)
    ! write(*,*) "wet_snow",wet_snow(448, 155)

    melt_snow = 0.
    where (wet_snow) melt_snow = (MELTws+max(0.,-PREFRws))
    where (dry_snow) melt_snow = (MELTds+max(0.,-PREFRds))
    where (new_snow) melt_snow = (MELTns+max(0.,-PREFRns))
    ! write(*,*) "melt_snow",melt_snow(448, 155)

    SNH_OUT(:,:,k)  = max(0.,tmpSNH + (snow)*30.5)
    ! write(*,*) "SNH_OUT",SNH_OUT(448, 155,k)
    ! write(*,*) "tmpSNH",tmpSNH(448, 155)
    ! write(*,*) "(snow(:,:,k))*30.5",(snow(448, 155))

    REFR_OUT(:,:,k)  = 0.
    where (wet_snow) REFR_OUT(:,:,k) = max(0.,PREFRws)
    where (dry_snow) REFR_OUT(:,:,k) = max(0.,PREFRds)
    where (new_snow) REFR_OUT(:,:,k) = max(0.,PREFRns)
    where (mask_in(:,:,k))  REFR_OUT(:,:,k) = min(.6*(SNH_OUT(:,:,k)/30.5), REFR_OUT(:,:,k))
    ! write(*,*) "REFR_OUT",REFR_OUT(448, 155,k)

    ! ice_melt_fract: to dertermine the percentage of ice melt in total melt
    ice_melt = (melt_snow - REFR_OUT(:,:,k) - tmpSNH/30.5 - snow)
    where ( ice_melt > 0.)
      ice_melt_fract  = ice_melt/melt_snow
    elsewhere
      ice_melt_fract  = 0.
    end where
    ! ice_melt_fract(:,:,k)  = (melt_snow(:,:,k)-REFR_OUT(:,:,k)-tmpSNH_in/30.5-snow(:,:,k))/melt_snow(:,:,k)     ! seems to work for matlab, but may involve a division by zero!!!
      ! write(*,*) "(melt_snow-REFR(:,:,month)-tmpSNH/30.5-snow)/melt_snow=",(melt_snow(448, 155,k)-REFR_OUT(448, 155,k)-tmpSNH_in(448, 155)/30.5-snow(448, 155,k))/melt_snow(448, 155,k)
      ! write(*,*) "(melt_snow-REFR(:,:,month)-tmpSNH/30.5-snow)",(melt_snow(448, 155,k)-REFR_OUT(448, 155,k)-tmpSNH_in(448, 155)/30.5-snow(448, 155,k))
      ! write(*,*) "(melt_snow-REFR(:,:,month))",(melt_snow(448, 155,k)-REFR_OUT(448, 155,k))
      ! write(*,*) "(tmpSNH/30.5)",(tmpSNH_in(448, 155)/30.5)
      ! write(*,*) "(snow(448, 155))",(snow(448, 155,k))
      ! write(*,*) "ice_melt_fract",ice_melt_fract(448, 155)
      ! write(*,*) "tmpSNH",tmpSNH(448, 155)
      ! write(*,*) "lastyear",lastyear(448, 155)

    MELT_OUT(:,:,k) = melt_snow*(1. - ice_melt_fract) + MELTi*ice_melt_fract
    ! write(*,*) "MELTi",MELTi(448, 155,k)
    ! write(*,*) "MELT_OUT",MELT_OUT(448, 155,k)

    SNH_OUT(:,:,k)  = max(0., tmpSNH + (snow - melt_snow + REFR_OUT(:,:,k))*30.5)
    ! write(*,*) "SNH_OUT",SNH_OUT(448, 155,k)

    ! Asnow(:,:,k)  = 0.
    where (wet_snow)
      Asnow = Aws
    ! write(*,*) "Ans",Asnow(448, 155,k)
    elsewhere (dry_snow)
      Asnow = Ads
    elsewhere (new_snow)
      Asnow = Ans
    elsewhere
      Asnow = 0.
    end where
    where (mask_in(:,:,k)) A_OUT(:,:,k)   = (Asnow * (1 - ice_melt_fract) + Ai * ice_melt_fract)
  !   write(*,*) "Asnow",Asnow(448, 155)
  ! write(*,*) "A_OUT",A_OUT(448, 155,k)

    where (mask_in(:,:,k)) ACC_OUT(:,:,k) = snow
    ! write(*,*) "ACC_OUT",ACC_OUT(448, 155,k)

    if (k==9) then
        tmpSNH = max(0.,SNH_OUT(:,:,k) - lastyear)
        lastyear = tmpSNH
    else
        tmpSNH = SNH_OUT(:,:,k)
    end if

    where (mask_in(:,:,k)) SMB_OUT(:,:,k) = (ACC_OUT(:,:,k) - MELT_OUT(:,:,k) + REFR_OUT(:,:,k))
    ! write(*,*) "SMB_OUT",SMB_OUT(448, 155,k)

  end do

  ! convert output data to expected units
  ! SNH: convert units from "mm" to "m"
  ! SMB: convert units from "mm day-1" to "kg m-2 second-1"
  ! ACC: convert units from "mm day-1" to "kg m-2 second-1"
  ! MELT: convert units from "mm day-1" to "kg m-2 second-1"
  ! REFR: convert units from "mm day-1" to "kg m-2 second-1"
  SNH_OUT  = SNH_OUT/1000.
  SMB_OUT  = SMB_OUT/(24.*60.*60.) 
  ACC_OUT  = ACC_OUT/(24.*60.*60.)
  MELT_OUT = MELT_OUT/(24.*60.*60.)
  REFR_OUT = REFR_OUT/(24.*60.*60.)

END SUBROUTINE dEBM_core


END MODULE MOD_MAIN
