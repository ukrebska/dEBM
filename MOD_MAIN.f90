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

  integer, intent(in) :: days_in
  real, intent(in) :: obl_in
  real, intent(out) :: DECL_out

  DECL_out = obl_in * sin(pi/180*360/365*(days_in-79))   ! 79 is the number of days since New Year for the spring equinox

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
! write(*,*) "elev", elev1(376,240)

! ! elev new
  ! elev1 = (decl + (90. - lat_in))*pi/180.
!   ! write(*,*) "elev1", elev1(376,240)
!   ! write(*,*) "elev_in", elev_in

  ! !products of sin(lat)*sin(decl) and cos(lat)*cos(decl)  is needed more than once
  sinphisind       = sin(pi/180*lat_in)*sin(pi/180*decl)
  cosphicosd       = cos(pi/180*lat_in)*cos(pi/180*decl)
  ! write(*,*) "sinphisind", sinphisind(376,240)
  ! write(*,*) "cosphicosd", cosphicosd(376,240)

  ! hourangle of sun rise and factor which relates solar flux density and
  ! daily insolation SW
   ! ha_0             = real(acos(-sinphisind/cosphicosd))
   ! sol_flux_fact_0  = (sinphisind*ha_0+cosphicosd*sin(ha_0))/pi
   ! write(*,*) "ha_0", ha_0(376,240)
   ! write(*,*) "sol_flux_fact_0", sol_flux_fact_0(376,240)

   ! ha_0 = 0.
   ! sol_flux_fact_0  = 0.
   ! where (abs(sinphisind/cosphicosd)<1.)
   ha_0 = -sinphisind/cosphicosd
     where (ha_0 > 1.) ha_0 = 1.
     where (ha_0 < -1.) ha_0 = -1.
   ha_0 = acos(ha_0)
  sol_flux_fact_0  = (sinphisind*ha_0+cosphicosd*sin(ha_0))/pi
   ! end where
   ! write(*,*) "-sinphisind/cosphicosd", -sinphisind(376,240)/cosphicosd(376,240)
   ! write(*,*) "ha_0", ha_0(376,240)
   ! write(*,*) "sol_flux_fact_0", sol_flux_fact_0(376,240)


  ! hourangle of sun rising above elev
  ha_elev = (sin(elev1)-sinphisind)/cosphicosd
    where (ha_elev > 1.) ha_elev = 1.
    where (ha_elev < -1.) elev1 = -1.
  ha_elev = acos(ha_elev)
    ! write(*,*) "(sin(elev1))", (sin(elev1))
    ! write(*,*) "(sin(elev1)-sinphisind)", (sin(elev1)-sinphisind(376,240))
    ! write(*,*) "(sin(elev1)-sinphisind)/cosphicosd", (sin(elev1)-sinphisind(376,240))/cosphicosd(376,240)
    ! write(*,*) "acos((sin(elev1)-sinphisind)/cosphicosd)", acos((sin(elev1)-sinphisind(376,240))/cosphicosd(376,240))
    ! ha_elev(:,:)           = acos(0.8)
    ! write(*,*) "ha_elev0", ha_elev0(376,240)
    ! write(*,*) "ha_elev", ha_elev(376,240)
    ! write(*,*) "ha_elev", maxval(ha_elev)
    ! write(*,*) "ha_elev", minval(ha_elev)


  ! duration of melt period in hours
  HOURS_OUT         = ha_elev*24/pi
  ! write(*,*) "HOURS_OUT", HOURS_OUT(376,240)
  ! write(*,*) "HOURS_OUT", minval(HOURS_OUT)

  ! proportion of daily insolation which is effective during melt period
  sol_flux_fact_e   = (sinphisind*ha_elev+cosphicosd*sin(ha_elev))/pi
  Q_OUT             = sol_flux_fact_e/sol_flux_fact_0

  ! where (sol_flux_fact_0==0.)
  !   Q_OUT = 0.
  ! elsewhere
  !   Q_OUT = sol_flux_fact_e/sol_flux_fact_0
  ! end where

  ! write(*,*) "sol_flux_fact_e", sol_flux_fact_e(376,240)
  ! write(*,*) "sol_flux_fact_0", sol_flux_fact_0(376,240)
  ! write(*,*) "Q_OUT", Q_OUT(376,240)

  FLUXFAC_OUT       = sol_flux_fact_0
  ! write(*,*) "FLUXFAC_OUT", FLUXFAC_OUT(376,240)

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



SUBROUTINE dEBMmodel_cloud(A_in, swd_in, Temp_in, PDD_in, rain_in, hours_in, q_in, MC_in, fluxfac_in, taucs_in, tauoc_in, MELT_OUT, REFR_OUT)
  ! ************************************************************************
  ! * calculates surface melt and refreeze                                 *
  ! ************************************************************************
  ! * INPUT: A,swd,T,PDD,rain,hours,q,MC,fluxfac                           *
  ! *   A             : Albedo                                             *
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
  ! ************************************************************************

  real, intent(in) :: A_in, taucs_in, tauoc_in
  real, intent(in), dimension(:,:)    :: swd_in, Temp_in, PDD_in, rain_in, hours_in, q_in, fluxfac_in
  logical, intent(in), dimension(:,:) :: MC_in
  real, intent(out), dimension(:,:)   :: MELT_OUT, REFR_OUT

  real, parameter :: Lf  = 3.34e5         ! Lf is latent heat of fusion
  real, parameter :: S0  = 1330                    ! solar constant
  real, allocatable, dimension(:,:)   :: Scs, Soc, Scsr, gwi
  real, allocatable, dimension(:,:)   :: meltcs, meltFDcs, meltFDoc, MELT_cs, MELT_oc

  allocate (Scs(xlen,ylen))
  allocate (Soc(xlen,ylen))
  allocate (Scsr(xlen,ylen))
  allocate (gwi(xlen,ylen))
  allocate (meltcs(xlen,ylen))
  allocate (meltFDcs(xlen,ylen))
  allocate (meltFDoc(xlen,ylen))
  allocate (MELT_cs(xlen,ylen))
  allocate (MELT_oc(xlen,ylen))

  Scs = fluxfac_in*S0*taucs_in        ! expected sunshine for clear sky conditions
  Soc = fluxfac_in*S0*tauoc_in        ! expected sunshine for completely overcast conditions
  ! write(*,*) "fluxfac_in", fluxfac_in(376,240)
  ! write(*,*) "Scs", Scs(376,240)
  ! write(*,*) "Soc", Soc(376,240)

  ! good weather index
  ! percentage of days in a month, which has clear sky conditions
  ! (there is just overcast or clear conditions, nothing inbetween)
  gwi = Scs*0
  where (Scs>50)
    gwi = min(1., (swd_in - Soc)/(Scs - Soc))
  end where
  ! write(*,*) "gwi", gwi(376,240)
!
  ! ?????
  ! why use Scsr, the maximum of real and expected sunshine, rather than real
  Scsr = max(Scs,swd_in)
  ! write(*,*) "swd_in", swd_in(376,240)
  ! write(*,*) "Scsr", Scsr(376,240)

  ! 0.05 is penetration
  meltcs = 0.
  where (MC_in) meltcs   = (((1-A_in)*Scsr*24*q_in)+(c2cs+c1cs*PDD_in)*hours_in)/Lf*60*60
  ! write(*,*) "A_in", A_in
  ! write(*,*) "Scsr", Scsr(376,240)
  ! write(*,*) "q_in", q_in(376,240)
  ! write(*,*) "PDD_in", PDD_in(376,240)
  ! write(*,*) "hours_in", hours_in(376,240)
  ! write(*,*) "MC_in", MC_in(376,240)
  ! write(*,*) "meltcs", meltcs(376,240)

  meltFDcs = 0.
  where (MC_in) meltFDcs = ((1-A_in)*Scsr+(c2cs+c1cs*Temp_in))/Lf*24*60*60
  meltFDoc = 0.
  where (MC_in) meltFDoc = ((1-A_in-.05)*Soc+(c2oc+(c1oc)*Temp_in))/Lf*24*60*60        ! we don't distinguish night and day here...
  ! MELT_cs
  where (isnan(meltcs))
    MELT_cs = 0.
  elsewhere
    MELT_cs = gwi*max(0.,meltcs)
  end where
  ! MELT_oc
  where (isnan(gwi*(meltFDcs-meltcs)+(1-gwi)*meltFDoc))
    MELT_oc = 0.
  elsewhere
    MELT_oc = max(0.,gwi*(meltFDcs-meltcs)+(1-gwi)*meltFDoc)
  end where
  !
  MELT_OUT = MELT_cs + MELT_oc
  ! write(*,*) "MELT_cs", MELT_cs(376,240)
  ! write(*,*) "MELT_oc", MELT_oc(376,240)

  ! write(*,*) "MC_in", MC_in(376,240)
  ! write(*,*) "meltcs", meltcs(376,240)
  ! write(*,*) "meltFDcs", meltFDcs(376,240)
  ! write(*,*) "meltFDoc", meltFDoc(376,240)
  ! write(*,*) "MELT_OUT", MELT_OUT(376,240)

  ! REFR_OUT = max(0.,min(MELT_OUT + rain_in, gwi*(meltcs-meltFDcs)))
  where (isnan(gwi*(meltcs-meltFDcs)))
    REFR_OUT = MELT_OUT + rain_in
  elsewhere (isnan(MELT_OUT + rain_in))
    REFR_OUT = gwi*(meltcs-meltFDcs)
  elsewhere
    REFR_OUT = min(MELT_OUT + rain_in, gwi*(meltcs-meltFDcs))
  end where
  !
  where (isnan(REFR_OUT))
    REFR_OUT = 0.
  elsewhere
    REFR_OUT = max(0.,REFR_OUT)
  end where

  ! write(*,*) "REFR_OUT", REFR_OUT(376,240)

END SUBROUTINE dEBMmodel_cloud




SUBROUTINE dEBM_core(tempm_in, radm_in, ppm_in, tmpSNH_in, lastyear_in, lat_in, mask_in, SNH_OUT, SMB_OUT, MELT_OUT, ACC_OUT, REFR_OUT, A_OUT)
  ! ************************************************************************
  ! * This SUBTROUTINE calculates the surface mass blance                  *
  ! ************************************************************************
  ! * INPUT: tempm, radm, ppm, tmpSNH, lastyear, lat, msk                  *
  ! *  tempm,radm,ppm are 3-dimension                                      *
  ! *  tmpSNH, lastyear,lat,msk are 2 dimension                            *
  ! *   tempm         : temperature                                        *
  ! *   radm          : daily short wave insolation at the surface         *
  ! *   ppm           : precipitation                                      *
  ! *   tmpSNH        : snow height of the last timestep                   *
  ! *   lastyear      : snow height at the end of last summer              *
  ! *   lat           : latitude                                           *
  ! *   mask          : mapping                                            *
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
  real, intent(out), dimension(:,:,:) :: SNH_OUT, SMB_OUT, MELT_OUT, ACC_OUT, REFR_OUT, A_OUT

  integer :: k
  real, allocatable, dimension(:,:,:) :: hoursns, qns, fluxfacns,&
                                          &hoursds, qds, fluxfacds,&
                                          &hoursws, qws, fluxfacws,&
                                          &hoursi, qi, fluxfaci,&
                                          &PDD0
  real, allocatable, dimension(:,:,:) :: solid0, solid, rain, snow
  logical, allocatable, dimension(:,:,:) :: tempmask
  real, allocatable, dimension(:,:,:) :: MELTns, PREFRns
  real, allocatable, dimension(:,:,:) :: MELTds, PREFRds
  real, allocatable, dimension(:,:,:) :: MELTws, PREFRws
  real, allocatable, dimension(:,:,:) :: MELTi, PREFRi
  logical, allocatable, dimension(:,:,:) :: new_snow, dry_snow, wet_snow
  real, allocatable, dimension(:,:,:) :: melt_snow, ice_melt_fract
  real, allocatable, dimension(:,:,:) :: Asnow
  real, allocatable, dimension(:,:)   :: tmpSNH, lastyear

  !
  allocate (hoursns(xlen,ylen,mlen))
  allocate (qns(xlen,ylen,mlen))
  allocate (fluxfacns(xlen,ylen,mlen))
  !
  allocate (hoursds(xlen,ylen,mlen))
  allocate (qds(xlen,ylen,mlen))
  allocate (fluxfacds(xlen,ylen,mlen))
  !
  allocate (hoursws(xlen,ylen,mlen))
  allocate (qws(xlen,ylen,mlen))
  allocate (fluxfacws(xlen,ylen,mlen))
  !
  allocate (hoursi(xlen,ylen,mlen))
  allocate (qi(xlen,ylen,mlen))
  allocate (fluxfaci(xlen,ylen,mlen))
  !
  allocate (PDD0(xlen,ylen,mlen))
  !
  allocate (solid0(xlen,ylen,mlen))
  allocate (solid(xlen,ylen,mlen))
  allocate (rain(xlen,ylen,mlen))
  allocate (snow(xlen,ylen,mlen))
  !
  allocate (tempmask(xlen,ylen,mlen))
  !
  allocate (MELTns(xlen,ylen,mlen))
  allocate (PREFRns(xlen,ylen,mlen))
  allocate (MELTds(xlen,ylen,mlen))
  allocate (PREFRds(xlen,ylen,mlen))
  allocate (MELTws(xlen,ylen,mlen))
  allocate (PREFRws(xlen,ylen,mlen))
  allocate (MELTi(xlen,ylen,mlen))
  allocate (PREFRi(xlen,ylen,mlen))
  !
  allocate (new_snow(xlen,ylen,mlen))
  allocate (dry_snow(xlen,ylen,mlen))
  allocate (wet_snow(xlen,ylen,mlen))
  allocate (melt_snow(xlen,ylen,mlen))
  !
  allocate (ice_melt_fract(xlen,ylen,mlen))
  allocate (Asnow(xlen,ylen,mlen))
  allocate (tmpSNH(xlen,ylen))
  allocate (lastyear(xlen,ylen))

  ! write(*,*) "calculating parameters c1 anc c2..."
  CALL CALC_PARAMETER(c1cs, c2cs, c1oc, c2oc)
  ! write(*,*) c1cs, c2cs, c1oc, c2oc

  ! write(*,*) "calculating angles...."
  CALL CALC_ANGLES(winkelns, winkelds, winkelws, winkeli)
  ! write(*,*) winkelns, winkelds, winkelws, winkeli

  !
  SMB_OUT    = tempm_in*0
  SNH_OUT    = tempm_in*0
  ! MELT_OUT   = tempm_in*0
  ! ACC_OUT    = tempm_in*0
  ! REFR_OUT   = tempm_in*0
  A_OUT      = tempm_in*0

  tmpSNH = tmpSNH_in
  lastyear = lastyear_in

  do k = 1, 12

    write(*,*) "calculates month",k

    ! write(*,*) "calculates the time that the sun is above a certain elevation angle...."
    CALL dEBM_sunny_hours_c(k, winkelns, lat_in(:,:,k), obliquity, hoursns(:,:,k), qns(:,:,k), fluxfacns(:,:,k))
    CALL dEBM_sunny_hours_c(k, winkelds, lat_in(:,:,k), obliquity, hoursds(:,:,k), qds(:,:,k), fluxfacds(:,:,k))
    CALL dEBM_sunny_hours_c(k, winkelws, lat_in(:,:,k), obliquity, hoursws(:,:,k), qws(:,:,k), fluxfacws(:,:,k))
    CALL dEBM_sunny_hours_c(k, winkeli, lat_in(:,:,k), obliquity, hoursi(:,:,k), qi(:,:,k), fluxfaci(:,:,k))

    ! write(*,*) k, " hoursns", hoursns(376,240,k)
    ! write(*,*) k, " qns ", qns(376,240,k)
    ! write(*,*) k, " fluxfacns ", fluxfacns(376,240,k)

    ! to determine the fraction of solid and liquid water in precipitation
    solid0(:,:,k) = 0.
    where (tempm_in(:,:,k)>-slim) solid0(:,:,k) = sin(tempm_in(:,:,k)/2/slim*pi)+1
    solid0(:,:,k) = 1-.5*solid0(:,:,k)
    where (tempm_in(:,:,k)<slim) solid(:,:,k) = solid0(:,:,k)
    solid(:,:,k) = max(0.,solid(:,:,k))

    ! write(*,*) "solid",solid(376,240,k)
    ! write(*,*) "ppm_in",ppm_in(376,240,k)

    snow(:,:,k) = 0.
    where (mask_in(:,:,k)) snow(:,:,k)  = solid(:,:,k) * ppm_in(:,:,k)            ! snow&rain  fractionation acc. to Pipes & Quick 1977
    rain(:,:,k) = 0.
    where (mask_in(:,:,k)) rain(:,:,k)  = (1 - solid(:,:,k)) * ppm_in(:,:,k)
    tempmask(:,:,k) = ((tempm_in(:,:,k) > -6.5) .AND. (mask_in(:,:,k)))
    ! write(*,*) "temp",tempm_in(376,240,k)
    ! write(*,*) "ppm",ppm_in(376,240,k)
    ! write(*,*) "radm_in",radm_in(376,240,k)
    ! write(*,*) "snow",snow(376,240,k)
    ! write(*,*) "rain",rain(376,240,k)
    ! write(*,*) "mask_in",mask_in(376,240,k)
    ! write(*,*) "tempmask",tempmask(376,240,k)

    ! PDD=PDD4(temp,3.5)*msk
    ! write(*,*) "calculates PDD...."
    CALL PDD4(tempm_in(:,:,k), stddev, PDD0(:,:,k))
    ! write(*,*) "PDD0",PDD0(376,240,k)

    ! write(*,*) "calculates dEBMmodel_cloud...."
    ! write(*,*) "Ans"
    CALL dEBMmodel_cloud(Ans, radm_in(:,:,k), tempm_in(:,:,k), PDD0(:,:,k), rain(:,:,k), hoursns(:,:,k), qns(:,:,k), tempmask(:,:,k), fluxfacns(:,:,k), taucs, tauoc, MELTns(:,:,k), PREFRns(:,:,k))
    ! write(*,*) "MELTns",MELTns(376,240,k)
    ! write(*,*) "PREFRns",PREFRns(376,240,k)
  ! write(*,*) "Ads"
    CALL dEBMmodel_cloud(Ads, radm_in(:,:,k), tempm_in(:,:,k), PDD0(:,:,k), rain(:,:,k), hoursds(:,:,k), qds(:,:,k), tempmask(:,:,k), fluxfacds(:,:,k), taucs, tauoc, MELTds(:,:,k), PREFRds(:,:,k))
    ! write(*,*) "Aws"
    CALL dEBMmodel_cloud(Aws, radm_in(:,:,k), tempm_in(:,:,k), PDD0(:,:,k), rain(:,:,k), hoursws(:,:,k), qws(:,:,k), tempmask(:,:,k), fluxfacws(:,:,k), taucs, tauoc, MELTws(:,:,k), PREFRws(:,:,k))
    ! write(*,*) "Ai"
    CALL dEBMmodel_cloud(Ai , radm_in(:,:,k), tempm_in(:,:,k), PDD0(:,:,k), rain(:,:,k), hoursi(:,:,k), qi(:,:,k), tempmask(:,:,k), fluxfaci(:,:,k), taucs, tauoc, MELTi(:,:,k), PREFRi(:,:,k))

    ! write(*,*) k, " hoursns", hoursns(376,240,k)
    ! write(*,*) k, " hoursds", hoursds(376,240,k)
    ! write(*,*) k, " hoursws", hoursws(376,240,k)
    ! write(*,*) k, " hoursi", hoursi(376,240,k)

    ! write(*,*) k, " fluxfacns", fluxfacns(376,240,k)
    ! write(*,*) k, " fluxfacds", fluxfacds(376,240,k)
    ! write(*,*) k, " fluxfacws", fluxfacws(376,240,k)
    ! write(*,*) k, " fluxfaci", fluxfaci(376,240,k)

    ! write(*,*) "calculates snow type...."
    new_snow(:,:,k) = (MELTns(:,:,k) < snow(:,:,k))
    dry_snow(:,:,k) = ((.NOT.new_snow(:,:,k)) .AND. (PREFRds(:,:,k) >= (MELTds(:,:,k) + rain(:,:,k))))
    wet_snow(:,:,k) = ((.NOT.new_snow(:,:,k)) .AND. (.NOT.dry_snow(:,:,k)))
    ! write(*,*) "MELTns",MELTns(376,240,k)
    ! write(*,*) "snow",snow(376,240,k)
  ! write(*,*) "new_snow",new_snow(376,240,k)
    ! write(*,*) "dry_snow",dry_snow(376,240,k)
    ! write(*,*) "wet_snow",wet_snow(376,240,k)

    melt_snow(:,:,k) = 0.
    where (wet_snow(:,:,k)) melt_snow(:,:,k) = (MELTws(:,:,k)+max(0.,-PREFRws(:,:,k)))
    where (dry_snow(:,:,k)) melt_snow(:,:,k) = (MELTds(:,:,k)+max(0.,-PREFRds(:,:,k)))
    where (new_snow(:,:,k)) melt_snow(:,:,k) = (MELTns(:,:,k)+max(0.,-PREFRns(:,:,k)))
    ! write(*,*) "melt_snow",melt_snow(:,:,k)

    REFR_OUT(:,:,k)  = 0.
    where (wet_snow(:,:,k)) REFR_OUT(:,:,k) = max(0.,PREFRws(:,:,k))
    where (dry_snow(:,:,k)) REFR_OUT(:,:,k) = max(0.,PREFRds(:,:,k))
    where (new_snow(:,:,k)) REFR_OUT(:,:,k) = max(0.,PREFRns(:,:,k))
    where (mask_in(:,:,k))  REFR_OUT(:,:,k) = min(.6*(tmpSNH/mth_days(k)+snow(:,:,k)), REFR_OUT(:,:,k))
    ! write(*,*) "REFR_OUT",REFR_OUT(376,240,k)

    ! ice_melt_fract: to dertermine the percentage of ice melt in total melt
    where ((melt_snow(:,:,k) - REFR_OUT(:,:,k) - tmpSNH/mth_days(k) - snow(:,:,k)) > 0.)
      ice_melt_fract(:,:,k)  = (melt_snow(:,:,k) - REFR_OUT(:,:,k) - tmpSNH/mth_days(k) - snow(:,:,k))/melt_snow(:,:,k)
    elsewhere
      ice_melt_fract(:,:,k)  = 0.
    end where
    ! ice_melt_fract(:,:,k)  = (melt_snow(:,:,k)-REFR_OUT(:,:,k)-tmpSNH_in/30.5-snow(:,:,k))/melt_snow(:,:,k)     ! seems to work for matlab, but may involve a division by zero!!!
      ! write(*,*) "(melt_snow-REFR(:,:,month)-tmpSNH/30.5-snow)/melt_snow=",(melt_snow(376,240,k)-REFR_OUT(376,240,k)-tmpSNH_in(376,240)/30.5-snow(376,240,k))/melt_snow(376,240,k)
      ! write(*,*) "(melt_snow-REFR(:,:,month)-tmpSNH/30.5-snow)",(melt_snow(376,240,k)-REFR_OUT(376,240,k)-tmpSNH_in(376,240)/30.5-snow(376,240,k))
      ! write(*,*) "(melt_snow-REFR(:,:,month))",(melt_snow(376,240,k)-REFR_OUT(376,240,k))
      ! write(*,*) "(tmpSNH/30.5)",(tmpSNH_in(376,240)/30.5)
      ! write(*,*) "(snow(376,240))",(snow(376,240,k))
      ! write(*,*) "ice_melt_fract",ice_melt_fract(576,491,k)

    MELT_OUT(:,:,k) = melt_snow(:,:,k)*(1 - ice_melt_fract(:,:,k)) + MELTi(:,:,k)*ice_melt_fract(:,:,k)
    ! write(*,*) "MELTi",MELTi(376,240,k)
    ! write(*,*) "MELT_OUT",MELT_OUT(376,240,k)

    SNH_OUT(:,:,k)  = max(0.,tmpSNH+(snow(:,:,k) - melt_snow(:,:,k) + REFR_OUT(:,:,k))*mth_days(k))
    ! write(*,*) "SNH_OUT",SNH_OUT(376,240,k)

    ! Asnow(:,:,k)    = -Ans*new_snow(:,:,k) - Ads*dry_snow(:,:,k) - Aws*wet_snow(:,:,k)
    ! write(*,*) "Asnow",Asnow(576,491,k)

    Asnow(:,:,k)  = 0.
    where (wet_snow(:,:,k))
      Asnow(:,:,k) = Aws
    ! write(*,*) "Ans",Asnow(376,240,k)
    elsewhere (dry_snow(:,:,k))
      Asnow(:,:,k) = Ads
    ! write(*,*) "Ans",Asnow(376,240,k)
    elsewhere (new_snow(:,:,k))
      Asnow(:,:,k) = Ans
    elsewhere
      Asnow(:,:,k) = 0.
    end where
    ! write(*,*) "Asnow",Asnow(576,491,k)

    where (mask_in(:,:,k)) A_OUT(:,:,k)   = (Asnow(:,:,k) * (1 - ice_melt_fract(:,:,k)) + Ai * ice_melt_fract(:,:,k))
    ! write(*,*) "A_OUT",A_OUT(576,491,k)

    where (mask_in(:,:,k)) ACC_OUT(:,:,k) = snow(:,:,k)
    ! write(*,*) "ACC_OUT",ACC_OUT(376,240,k)

    if (k==9) then
        tmpSNH = max(0.,SNH_OUT(:,:,k) - lastyear)
        lastyear = tmpSNH
    else
        tmpSNH = SNH_OUT(:,:,k)
    end if

    where (mask_in(:,:,k)) SMB_OUT(:,:,k) = (ACC_OUT(:,:,k) - MELT_OUT(:,:,k) + REFR_OUT(:,:,k))

  end do

END SUBROUTINE dEBM_core


END MODULE MOD_MAIN
