MODULE MOD_MAIN
  ! ************************************************************************
  ! * This MODULE calculates surface mass balance                          *
  ! ************************************************************************
  ! * MOD_MAIN                                                             *
  ! *     | dEBM_core                                                      *
  ! *         | dEBM_decl                                                  *
  ! *         | dEBM_fluxfac                                               *
  ! *         | dEBM_sunny_hours_c                                         *
  ! *         | PDD4                                                       *
  ! *         | dEBMmodel_fullrad                                          *
  ! ************************************************************************

  USE MOD_PRE
  USE MOD_DATA

  implicit none

  real(kind=WP), parameter :: pi     = 3.141592653
  real(kind=WP), dimension(12) :: mth = (/1,2,3,4,5,6,7,8,9,10,11,12/)
  real(kind=WP), dimension(12) :: mth_len = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  real(kind=WP), dimension(12) :: days = (/16,45,75,105,136,166,197,228,258,289,319,350/)

contains

SUBROUTINE dEBM_decl(days, obl, DECL)
  ! ************************************************************************
  ! * dEBM_decl calcualtes declination                                     *
  ! *                                                                      *
  ! *     assumes eccentricity is 1                                        *
  ! *     spring equinox is fixed on March 20th                            *
  ! *     neglects leap years                                              *
  ! ************************************************************************
  ! * INPUT:                                                               *
  ! *   days           :  days in each month                               *
  ! *   obliqu         :  obliquity                                        *
  ! ************************************************************************
  ! * OUTPUT:                                                              *
  ! *   DECL          : declination                                        *
  ! ************************************************************************

  real(kind=WP), intent(in) :: days
  real(kind=WP), intent(in) :: obl
  real(kind=WP), intent(out) :: DECL

  DECL = obl * sin(pi/180.0_WP*360.0_WP/365.0_WP*(days-79.0_WP))   ! 79 is the number of days since New Year for the spring equinox

  ! debug
  if (debug_switch) then
    write(*,*) "surviving dEBM_decl"
  end if

END SUBROUTINE dEBM_decl


SUBROUTINE dEBM_fluxfac(mth, latm, obl, sol_flux_fact_0)
  ! ************************************************************************
  ! dEBM_fluxfac                                                           *
  !       calculates daily insolation                                      *
  ! ************************************************************************
  ! * INPUT:                                                               *
  ! *   mth           : months                                             *
  ! *   latm          : latitude                                           *
  ! *   obl           : obliquity                                          *
  ! ************************************************************************
  ! * OUTPUT:                                                              *
  ! *   sol_flux_fact_0  : daily insolation                                *
  ! ************************************************************************

  real(kind=WP), intent(in) :: obl
  real(kind=WP), intent(in), dimension(:,:) :: latm
  integer, intent(in) :: mth

  real(kind=WP), intent(out), dimension(:,:) :: sol_flux_fact_0

  real(kind=WP), allocatable, dimension(:,:) :: decl
  real(kind=WP), allocatable, dimension(:,:) :: sinphisind,cosphicosd
  real(kind=WP), allocatable, dimension(:,:) :: ha_0

  ! decl
  allocate (decl(xlen, ylen))
  CALL dEBM_decl(days(mth), obl, decl(1,1))
  do i = 1, xlen
      do j = 1, ylen
          decl(i, j) = decl(1,1)
      end do
  end do

  ! hourangle of sun rise
  allocate (sinphisind(xlen, ylen))
  allocate (cosphicosd(xlen, ylen))
  sinphisind = sin(pi/180.0_WP*latm)*sin(pi/180.0_WP*decl)
  cosphicosd = cos(pi/180.0_WP*latm)*cos(pi/180.0_WP*decl)
  allocate (ha_0(xlen, ylen))
   where (-sinphisind/cosphicosd > 1.0_WP)
     ha_0 = 0.0_WP  ! acos(1.)
   elsewhere (-sinphisind/cosphicosd < -1.0_WP)
     ha_0 = pi      ! acos(-1.)
   elsewhere
     ha_0 = acos(-sinphisind/cosphicosd)
   end where

   ! daily insolation SW
   sol_flux_fact_0  = (sinphisind*ha_0+cosphicosd*sin(ha_0))/pi

   ! deallocate
   deallocate(decl, sinphisind, cosphicosd, ha_0)

   ! debug
   if (debug_switch) then
     write(*,*) "surviving dEBM_fluxfac"
   end if

END SUBROUTINE dEBM_fluxfac


SUBROUTINE dEBM_sunny_hours_c(mth, elev, latm, obl, HOURS, Q, FLUXFAC)
  ! ************************************************************************
  ! dEBM_sunny_hours_c                                                     *
  !   calculates the hours that the sun is above a certain elevation angle *
  !   daily insolation and solar density                                   *
  ! ************************************************************************
  ! * INPUT:                                                               *
  ! *   mth      : month                                                   *
  ! *   elev     : elevation angle                                         *
  ! *   latm     : latitude                                                *
  ! *   obl      : obliquity                                               *
  ! ************************************************************************
  ! * OUTPUT:                                                              *
  ! *   HOURS    : duration of melt period in hours                        *
  ! *   Q        : ratio between insolation of melt period and whole day   *
  ! *   FLUXFAC  : daily insolation                                        *
  ! ************************************************************************

  real(kind=WP), intent(in) :: obl
  real(kind=WP), intent(in), dimension(:,:) :: elev, latm
  integer, intent(in) :: mth

  real(kind=WP), intent(out), dimension(:,:) :: HOURS, Q, FLUXFAC

  real(kind=WP), allocatable, dimension(:,:) :: elev1, decl
  real(kind=WP), allocatable, dimension(:,:) :: sinphisind,cosphicosd
  real(kind=WP), allocatable, dimension(:,:) :: ha_0,ha_elev0,ha_elev
  real(kind=WP), allocatable, dimension(:,:) :: sol_flux_fact_0,sol_flux_fact_e

  ! calc mid-month declination
  allocate (decl(xlen, ylen))
  CALL dEBM_decl(days(mth), obl, decl(1,1))
  ! inflating declination to ( xlen x ylen  )
  do i = 1, xlen
      do j = 1, ylen
          decl(i, j) = decl(1,1)
      end do
  end do

  ! elevation angle
  allocate (elev1(xlen, ylen))
  elev1(:,:) = elev*pi/180.0_WP

  ! hourangle of sun rise
  allocate (sinphisind(xlen, ylen))
  allocate (cosphicosd(xlen, ylen))
  sinphisind = sin(pi/180.0_WP*latm)*sin(pi/180.0_WP*decl)
  cosphicosd = cos(pi/180.0_WP*latm)*cos(pi/180.0_WP*decl)
  allocate (ha_0(xlen, ylen))
  allocate (ha_elev0(xlen, ylen))
  ha_elev0 = -sinphisind/cosphicosd
  where (ha_elev0 > 1.0_WP)
    ha_0 = 0.0_WP  ! acos(1.)
  elsewhere (ha_elev0 < -1.0_WP)
    ha_0 = pi      ! acos(-1.)
  elsewhere
    ha_0 = acos(ha_elev0)
  end where

  ! daily insolation
  allocate (sol_flux_fact_0(xlen, ylen))
  sol_flux_fact_0  = (sinphisind*ha_0+cosphicosd*sin(ha_0))/pi
  FLUXFAC = sol_flux_fact_0

  ! hourangle of sun rising above elev
  allocate (ha_elev(xlen, ylen))
  ha_elev0 = (sin(elev1)-sinphisind)/cosphicosd
  where (ha_elev0 > 1.0_WP)
    ha_elev = 0.0_WP  ! acos(1.)
  elsewhere (ha_elev0 < -1.0_WP)
    ha_elev = pi      ! acos(-1.)
  elsewhere
    ha_elev = acos(ha_elev0)
  endwhere

  ! duration of melt period
  HOURS = ha_elev*24.0_WP/pi

  ! proportion of daily insolation which is effective during melt period
  allocate (sol_flux_fact_e(xlen, ylen))
  sol_flux_fact_e = (sinphisind*ha_elev+cosphicosd*sin(ha_elev))/pi

  ! factor relates to solar flux density
  where (sol_flux_fact_0==0.0_WP)
    Q = 0.0_WP
  elsewhere
    Q = sol_flux_fact_e/sol_flux_fact_0
  end where

  ! deallocate
  deallocate(decl, elev1, sinphisind, cosphicosd)
  deallocate(ha_0, ha_elev0, ha_elev, sol_flux_fact_0,sol_flux_fact_e)

  ! debug
  if (debug_switch) then
    write(*,*) "surviving dEBM_sunny_hours_c"
  end if

END SUBROUTINE dEBM_sunny_hours_c



SUBROUTINE PDD4(temp, stddev, PDD)
  ! ************************************************************************
  ! * PDD is approximated as described in                                  *
  ! * Calov R and Greve R (2005)                                           *
  ! * Correspondence. A semi-analytical solution for the positive          *
  ! * degree-day model with stochastic temperature variations.             *
  ! * J. Glaciol., 51 (172), 173?175                                       *
  ! * (doi:10.3189/172756505781829601)                                     *
  ! *                                                                      *
  ! * INPUT:                                                               *
  ! *   temp          : surface temperature                                *
  ! *   stddev        : constant standard deveation                        *
  ! ************************************************************************
  ! * OUTPUT:                                                              *
  ! *   PDD           :                                                    *
  ! ************************************************************************

  real(kind=WP), intent(in) :: stddev
  real(kind=WP), intent(in), dimension(:,:) :: temp
  real(kind=WP), intent(out), dimension(:,:) :: PDD

  PDD = stddev/sqrt(2.0_WP*pi)*exp(-(temp**2.0_WP)/(2.0_WP*(stddev**2)))+temp/2.0_WP*erfc(-temp/(sqrt(2.0_WP)*stddev))

  ! debug
  if (debug_switch) then
    write(*,*) "surviving PDD4"
  end if

END SUBROUTINE PDD4



SUBROUTINE dEBMmodel_fullrad(A, Aoc, swd_cs, swd_oc, temp, pdd, rain, hours, q, c1cs, c2cs, c1oc, c2oc, MC, cc, MELT, REFR)
  ! ************************************************************************
  ! * dEBMmodel_fullrad calculates surface melt and refreeze               *
  ! ************************************************************************
  ! * INPUT:                                                               *
  ! *   A             : Albedo                                             *
  ! *   Aoc           : Albedo of overcast days                          *
  ! *   swd_cs        : surface short wave insolation for clearsky         *
  ! *   swd_oc        : surface short wave insolation for overcast       *
  ! *   temp          : surface temperature                                *
  ! *   pdd           : an approx. of melt-period temperature              *
  ! *   rain          : liquid precipitation                               *
  ! *   hours         : duration of melt period in hours                   *
  ! *   q             : ratio of effective radiation                       *
  ! *   MC            : the background melt condition (e.g. T>-6.5)        *
  ! *   cc            :  cloud cover                                       *
! **************************************************************************
  ! * OUTPUT:                                                              *
  ! *   MELT          : melt                                               *
  ! *   REFR          : refreeze                                           *
  ! ************************************************************************

  real(kind=WP), intent(in), dimension(:,:) :: A, Aoc
  real(kind=WP), intent(in), dimension(:,:) :: c1cs, c2cs, c1oc, c2oc
  real(kind=WP), intent(in), dimension(:,:)    :: swd_cs, swd_oc
  real(kind=WP), intent(in), dimension(:,:)    :: temp, pdd, rain, hours, q, cc
  logical, intent(in), dimension(:,:) :: MC

  real(kind=WP), intent(out), dimension(:,:)   :: MELT, REFR

  real(kind=WP), parameter :: Lf  = 3.34e5         ! Lf is latent heat of fusion
  real(kind=WP), allocatable, dimension(:,:)   :: meltcs, meltoc, refrcs, refroc, meltFDcs, meltFDoc

  allocate (meltcs(xlen, ylen))
  allocate (meltoc(xlen, ylen))
  allocate (refrcs(xlen, ylen))
  allocate (refroc(xlen, ylen))
  allocate (meltFDcs(xlen, ylen))
  allocate (meltFDoc(xlen, ylen))

  where (MC) meltcs   = (((1.0_WP-A)*swd_cs*24.0_WP*q) + (c2cs+c1cs*pdd)*hours)/Lf*60.0_WP*60.0_WP ! melt clearsky
  where (MC) meltFDcs = ((1.0_WP-A)*swd_cs + (c2cs+c1cs*temp))/Lf*24.0_WP*60.0_WP*60.0_WP          ! melt clearsky fullday
  where (MC) meltFDoc = ((1.0_WP-Aoc)*swd_oc + (c2oc+(c1oc)*temp))/Lf*24.0_WP*60.0_WP*60.0_WP      ! melt overcast fullday
                                                                                                   ! we don't distinguish night and day here...

  refroc = max(-meltFDoc, 0.0_WP)               ! refreeze overcast
  meltoc = max( meltFDoc, 0.0_WP)               ! melt overcast
  meltcs = max(0.0_WP, max(meltFDcs, meltcs))   ! melt clearsky
  refrcs = max(0.0_WP, meltcs - meltFDcs)       ! refreeze clearsky

  MELT = (1.0_WP-cc)*meltcs + cc*meltoc;
  REFR = min( MELT + rain, (1.0_WP-cc)*refrcs + cc*refroc )

  ! deallocate
  deallocate( meltcs, meltoc, refrcs, refroc, meltFDcs, meltFDoc )

  ! debug
  if (debug_switch) then
    write(*,*) "surviving dEBMmodel_fullrad"
  end if

END SUBROUTINE



SUBROUTINE dEBM_core(tempm, swdm, swd_TOAm, emiss, clcov, ppm, tmpSNH, lastyear, latm, mask, obl, mth_str, SNH, SMB, MELT, REFR, A, SNOW, RAIN, S07)
  ! ************************************************************************
  ! * dEBM_core calculates the surface mass blance                         *
  ! ************************************************************************
  ! input:                                                                 *
  ! tempm    : 2m temperature as 12 monthly means (downscaled)
  ! swdm     : downward shortwave radiation at surface (12 monthly means)
  ! swd_TOAm : downward shortwave radiation at top of atm. (12 monthly means)
  ! emiss    : atm. emissivity calculated on coarse grid,
  !            then interpolated (12 monthly means)
  ! clcov    : cloud cover [0 1]  (12 monthly means)
  ! ppm      : precipitation in mm/day (12  monthly means)
  ! tmpSNH  : last month's snow height (mm)
  ! lastyear: last September's snow height (mm)
  ! latm     : latitude
  ! mask     : ice mask
  ! obl      : obliquity
  ! mth_str  : list of month (first year [9:12 1:12], then [1:12])
  !*****************************************************************************
  ! output:
  ! SNH      : snow height in mm
  ! SMB      : surface mass balance (mm/day)
  ! MELT     : melt (mm/day)
  ! REFR     : refreezing (mm/day)
  ! A        : Albedo
  ! SNOW     : solid precipitation (mm/day)
  ! RAIN     : liquid precipitation (mm/day)
  ! S07      : summer solar density flux (TOA), just to check
  ! tmpSNH   : December month's snow height (mm)
  ! lastyear :  September's snow height (mm)
  !******************************************************************************

  integer, intent(in) :: mth_str
  real(kind=WP), intent(in) :: obl
  logical, intent(in), dimension(:,:,:) :: mask
  real(kind=WP), intent(in), dimension(:,:,:)    :: tempm, swdm, &
                                                      &emiss, clcov, ppm, latm

  real(kind=WP), intent(inout), dimension(:,:,:) :: swd_TOAm   ! if the swd_TOAm is not input
                                                               ! the input field will be zero, then recalculate
  real(kind=WP), intent(inout), dimension(:,:)   :: tmpSNH, lastyear ! snow height of December and September need to be saved
  real(kind=WP), intent(out), dimension(:,:,:)   :: SNH, SMB, MELT, REFR, &
                                                      &A, SNOW, RAIN
  real(kind=WP) , intent(out):: S07

  !  parameters
  ! TODO(sxu): better to be put in namelists
  real(kind=WP), parameter :: slim   = 7.       ! melt/refreeze threshold
  real(kind=WP), parameter :: Ansc    = .85     ! Albedo for new snow
  real(kind=WP), parameter :: Adsc    = .74      ! Albedo for dry snow
  real(kind=WP), parameter :: Awsc   = .65     ! Albedo for wet snow
  real(kind=WP), parameter :: Ai     = .5     ! Albedo for ice
  real(kind=WP), parameter :: sn_aging = .03  ! Albedo of wet snow darkens with temperature...
  real(kind=WP), parameter :: lhf    = 2;       ! residual heat flux (e.g. due to flux into the ice sheet
  real(kind=WP), parameter :: tau_cs = .8;      ! expected clear sky transmissivity
  ! emissivity depends on cloud cover and greenhouse gases (incl. water vapor)
  real(kind=WP), parameter :: epsa_cs = .72     ! eps_cs = emiss_cs + emiss_gas; (Sedlar & Hock, 2009)
  real(kind=WP), parameter :: epsa_oc = epsa_cs+.18     ! eps_oc = eps_cs + const;       (Konig-Langlo, 1994)
  real(kind=WP), parameter :: epsi   = .95     ! emissivity of ice
  real(kind=WP), parameter :: beta   = 10.     ! turbulent heat transfer coeff
  real(kind=WP), parameter :: bolz   = 5.67e-8 ! Stefan-Boltzmann constant
  real(kind=WP), parameter :: T0     = 273.15  ! melt point in K
  real(kind=WP), parameter :: Tmin   = -6.5    ! background melt condition

  integer :: month
  real :: ml
  real(kind=WP), allocatable, dimension(:,:) :: hoursns, qns, fluxfacns,&
                                          &hoursds, qds, fluxfacds,&
                                          &hoursws, qws, fluxfacws,&
                                          &temp, PDD
  real(kind=WP), allocatable, dimension(:,:) :: solid, cc,&
                                          &Ans, Ads, Aws,&
                                          &swd, swd_TOA, swd_cs, swd_oc
  real(kind=WP), allocatable, dimension(:,:) :: tau, tau_oc
  real(kind=WP), allocatable, dimension(:,:) :: epsa_csr, emiss_gas, epsa_cst, epsa_oct
  logical, allocatable, dimension(:,:) :: tmpmask
  real(kind=WP), allocatable, dimension(:,:) :: winkelns, winkelds, winkelws
  real(kind=WP), allocatable, dimension(:,:) :: c1cs, c2cs, c1oc, c2oc
  real(kind=WP), allocatable, dimension(:,:) :: MELTns, PREFRns
  real(kind=WP), allocatable, dimension(:,:) :: MELTds, PREFRds
  real(kind=WP), allocatable, dimension(:,:) :: MELTws, PREFRws, snh_est
  logical, allocatable, dimension(:,:) :: old_wet, new_snow, dry_snow, wet_snow
  real(kind=WP), allocatable, dimension(:,:,:) :: ff
  real(kind=WP), dimension(12) :: S0
  integer, dimension(3) :: min_lat_idx
  real :: swd_lat_min, fluxfactor

  allocate (hoursns(xlen, ylen))
  allocate (qns(xlen, ylen))
  allocate (fluxfacns(xlen, ylen))
  allocate (hoursds(xlen, ylen))
  allocate (qds(xlen, ylen))
  allocate (fluxfacds(xlen, ylen))
  allocate (hoursws(xlen, ylen))
  allocate (qws(xlen, ylen))
  allocate (fluxfacws(xlen, ylen))
  !
  allocate (temp(xlen, ylen))
  allocate (PDD(xlen, ylen))
  allocate (tmpmask(xlen, ylen))
  allocate (solid(xlen, ylen))
  !
  allocate (cc(xlen, ylen))
  allocate (Ans(xlen, ylen))
  allocate (Ads(xlen, ylen))
  allocate (Aws(xlen, ylen))

  allocate (swd(xlen, ylen))
  allocate (swd_TOA(xlen, ylen))
  allocate (swd_cs(xlen, ylen))
  allocate (swd_oc(xlen, ylen))

  allocate (tau(xlen, ylen))
  allocate (tau_oc(xlen, ylen))

  allocate (epsa_csr(xlen, ylen))
  allocate (emiss_gas(xlen, ylen))
  allocate (epsa_cst(xlen, ylen))
  allocate (epsa_oct(xlen, ylen))
  !
  allocate (MELTns(xlen, ylen))
  allocate (PREFRns(xlen, ylen))
  allocate (MELTds(xlen, ylen))
  allocate (PREFRds(xlen, ylen))
  allocate (MELTws(xlen, ylen))
  allocate (PREFRws(xlen, ylen))
  !
  allocate (old_wet(xlen, ylen))
  allocate (new_snow(xlen, ylen))
  allocate (dry_snow(xlen, ylen))
  allocate (wet_snow(xlen, ylen))
  !
  allocate (c1cs(xlen, ylen))
  allocate (c2cs(xlen, ylen))
  allocate (c1oc(xlen, ylen))
  allocate (c2oc(xlen, ylen))
  !
  allocate (snh_est(xlen, ylen))

  !Summer solar flux density
  S0(1:12)=1330.0_WP
  ! TODO:
  ! S0 should be a variable funct. f(eccentricity, longitude of perihelion,
  ! declination) for paleo applications we can calculate it from SWD_TOA
  ! ideally we should improve the calculations of the declination and S0

  allocate (ff(xlen, ylen, mlen))
  ff = 0.0_WP
  do month=1, 12
    CALL dEBM_fluxfac(month, latm(:,:,month), obl, ff(:,:,month))
  end do

  ! if swd_TOA is not given, S0 is constant 1330
  ! if swd_TOA is given, we estimate summer S0 from swd_TOA at 65N
  ! Here, we use 65N, because insolation at 65N can be compared to literature just to doublecheck
  ! TODO: at some point we should calculate this from orbital parameters
  if (use_shortwave_radiation_TOA) then
    min_lat_idx=minloc(abs(latm(:,:,:)-65.0_WP))     ! locate 65N
    do month=5, 9
      fluxfactor = ff(min_lat_idx(1),min_lat_idx(2),month)
      swd_lat_min = swd_TOAm(min_lat_idx(1),min_lat_idx(2),month);
      if (swd_lat_min>100.0_WP) then
        S0(month) = swd_lat_min/fluxfactor
      end if
    end do
  else
    do month = 1, 12
      swd_TOAm(:,:,month) = (S0(month)*ff(:,:,month))
    end do
  end if

  ! debug
  if (debug_switch) then
    write(*,*) "swd_TOAm",swd_TOAm(debug_lon, debug_lat,month)
  end if

  S07 = S0(7)
  ! initialization
  SMB(:,:,:)    = 0.0_WP
  SNH(:,:,:)    = 0.0_WP
  MELT(:,:,:)   = 0.0_WP
  REFR(:,:,:)   = 0.0_WP
  A(:,:,:)      = 0.0_WP
  RAIN(:,:,:)   = 0.0_WP
  SNOW(:,:,:)   = 0.0_WP
  wet_snow(:,:) = .FALSE.

  do month = mth_str, 12

    write(*,*) "calculates month", month

    ! TODO: set 30.5 for comparsion
    ! otherwise, use the real length of each month
    ! AND also make sure type of ml is integer! (Line 365-366)
    ! ml      = mth_len(month)
    ml      = 30.5
    cc      = clcov(:,:,month)
    swd     = swdm(:,:,month)
    swd_TOA = swd_TOAm(:,:,month)
    tmpmask = ((tempm(:,:,month) > -6.5) .AND. (mask(:,:,month)))
    CALL PDD4(temp, stddev, PDD)
    where (.NOT. tmpmask) PDD = 0.

    ! Albedo
    Ans = Ans
    Ads = Ads
    Aws = max(Ai, Awsc - sn_aging*PDD)
    ! debug
    if (debug_switch) then
      write(*,*) "cc",cc(debug_lon, debug_lat)
      write(*,*) "swd",swd(debug_lon, debug_lat)
      write(*,*) "swd_TOAm",swd_TOAm(debug_lon, debug_lat,month)
      write(*,*) "tempm",tempm(debug_lon, debug_lat,month)
      write(*,*) "mask",mask(debug_lon, debug_lat,month)
      write(*,*) "tmpmask",tmpmask(debug_lon, debug_lat)
      write(*,*) "PDD",PDD(debug_lon, debug_lat)
      write(*,*) "Aws",Aws(debug_lon, debug_lat)
    end if

    ! to determine the fraction of solid and liquid water in precipitation
    ! snow&rain  fractionation acc. to Pipes & Quick 1977
    ! solid
    where (tempm(:,:,month)>slim)
      solid = 0.0_WP
    elsewhere (tempm(:,:,month)<-slim)
      solid = 1.0_WP
    elsewhere
      solid = 1.0_WP - 0.5_WP*(sin(tempm(:,:,month)/2.0_WP/slim*pi) + 1.0_WP)
    end where
    where (mask(:,:,month)) SNOW(:,:,month)  = solid * ppm(:,:,month)
    where (mask(:,:,month)) RAIN(:,:,month)  = (1.0_WP - solid) * ppm(:,:,month)
    ! debug
    if (debug_switch) then
      write(*,*) "solid",solid(debug_lon, debug_lat)
      write(*,*) "SNOW",SNOW(debug_lon, debug_lat, month)
      write(*,*) "RAIN",RAIN(debug_lon, debug_lat, month)
    end if

    ! monthly mean transmissiv.
    ! in winter & avoid numeric problems when TOA ~ 0
    where (swd_TOA>5.0_WP)
      tau = min(0.9_WP, max( 0.1_WP, swd/swd_TOA ))
    elsewhere
      tau = 0.5_WP
    end where

    ! shortwave radiation of clear sky and overcast days
    ! on higher elevation it may happen that tau>tau_cs
    ! also we take care of extremely clear months (probably never a problem)
    swd_cs = max(tau_cs*S0(month)*ff(:,:,month), swd)
    where (cc>.05_WP)
      swd_oc = (swd - swd_cs*(1.0_WP-cc))/cc
    elsewhere
      swd_oc = swd_cs
    end where

    where (ff(:,:,month) /= 0.) tau_oc = min(tau_cs, swd_oc/(ff(:,:,month)*S0(month)))
    tau_oc   = min(0.7_WP, tau_oc)
    tau_oc   = max(0.0_WP, tau_oc)
    !
    if (debug_switch) then
      write(*,*) "swd_oc",swd_oc(debug_lon, debug_lat)
    end if

    ! clear sky and overcast emissivity
    epsa_csr  = epsa_cs ! min(epsa_cs,emiss(:,:,month))
    emiss_gas = emiss(:,:,month) - (1.0_WP-cc)*epsa_csr - cc*epsa_oc
    epsa_cst  = epsa_csr + emiss_gas
    epsa_oct  = epsa_oc  + emiss_gas
    !
    if (debug_switch) then
      write(*,*) "tau_oc",tau_oc(debug_lon, debug_lat)
      write(*,*) "epsa_oct",epsa_oct(debug_lon, debug_lat)
    end if

    ! c1, c2 are determined locally and monthly
    c2cs = (-epsi+epsa_cst*epsi)*bolz*(T0**4)-lhf
    c1cs = (epsi*epsa_cst*4.*bolz*(T0**3)+beta)
    c2oc = (-epsi+epsa_oct*epsi)*bolz*(T0**4)-lhf
    c1oc = (epsi*epsa_oct*4.*bolz*(T0**3)+beta)

    winkelns = asin(-c2cs/(1.0_WP-Ans)/(S0(month)*tau_cs))*180.0_WP/pi
    winkelds = asin(-c2cs/(1.0_WP-Ads)/(S0(month)*tau_cs))*180.0_WP/pi
    winkelws = asin(-c2cs/(1.0_WP-Aws)/(S0(month)*tau_cs))*180.0_WP/pi

    CALL dEBM_sunny_hours_c(month, winkelns, latm(:,:,month), obliquity, hoursns, qns, fluxfacns)
    CALL dEBM_sunny_hours_c(month, winkelds, latm(:,:,month), obliquity, hoursds, qds, fluxfacds)
    CALL dEBM_sunny_hours_c(month, winkelws, latm(:,:,month), obliquity, hoursws, qws, fluxfacws)
    !
    if (debug_switch) then
      write(*,*) "c2cs",c2cs(debug_lon, debug_lat)
      write(*,*) "c1cs",c1cs(debug_lon, debug_lat)
      write(*,*) "winkelns",winkelns(debug_lon, debug_lat)
      write(*,*) "winkelds",winkelds(debug_lon, debug_lat)
      write(*,*) "winkelws",winkelws(debug_lon, debug_lat)
      write(*,*) "hoursns",hoursns(debug_lon, debug_lat)
      write(*,*) "qns",qns(debug_lon, debug_lat)
      write(*,*) "fluxfacns",fluxfacns(debug_lon, debug_lat)
    end if

    ! MELT + REFR
    CALL dEBMmodel_fullrad(Ans, Ans+0.05_WP, swd_cs, swd_oc, tempm(:,:,month), PDD, RAIN(:,:,month), hoursns, qns, c1cs, c2cs, c1oc, c2oc, tmpmask, cc, MELTns, PREFRns)
    CALL dEBMmodel_fullrad(Ads, Ads+0.05_WP, swd_cs, swd_oc, tempm(:,:,month), PDD, RAIN(:,:,month), hoursds, qds, c1cs, c2cs, c1oc, c2oc, tmpmask, cc, MELTds, PREFRds)
    CALL dEBMmodel_fullrad(Aws, Aws+0.05_WP, swd_cs, swd_oc, tempm(:,:,month), PDD, RAIN(:,:,month), hoursws, qws, c1cs, c2cs, c1oc, c2oc, tmpmask, cc, MELTws, PREFRws)
    !
    if (debug_switch) then
      write(*,*) "RAIN",RAIN(debug_lon, debug_lat,month)
      write(*,*) "MELTns",MELTns(debug_lon, debug_lat)
      write(*,*) "PREFRns",PREFRns(debug_lon, debug_lat)
    end if

    !  snow type
    old_wet  = (wet_snow .AND. (PREFRws < (MELTws + RAIN(:,:,month))))
    new_snow = (MELTns <= SNOW(:,:,month))
    dry_snow = ((.NOT.new_snow) .AND. (PREFRds >= (MELTds + RAIN(:,:,month))) .AND. (.NOT.old_wet))
    wet_snow = ((.NOT.new_snow) .AND. ((.NOT.dry_snow) .OR. (old_wet)))
    !
    if (debug_switch) then
      write(*,*) "wet_snow",wet_snow(debug_lon, debug_lat)
      write(*,*) "dry_snow",dry_snow(debug_lon, debug_lat)
      write(*,*) "new_snow",new_snow(debug_lon, debug_lat)
    end if

    ! MELT
    where (wet_snow)   MELT(:,:,month) = (MELTws + max(0.0_WP,-PREFRws))
    where (dry_snow)   MELT(:,:,month) = (MELTds + max(0.0_WP,-PREFRds))
    where (new_snow)   MELT(:,:,month) = (MELTns + max(0.0_WP,-PREFRns))
    !
    if (debug_switch) then
      write(*,*) "MELTws",MELTws(debug_lon, debug_lat)
      write(*,*) "MELTds",MELTds(debug_lon, debug_lat)
      write(*,*) "MELTns",MELTns(debug_lon, debug_lat)
      write(*,*) "MELT",MELT(debug_lon, debug_lat,month)
    end if

    ! REFR
    snh_est = max(0.0_WP, tmpSNH + SNOW(:,:,month)*ml)
    where (wet_snow) REFR(:,:,month) = max(0.0_WP,PREFRws)
    where (dry_snow) REFR(:,:,month) = max(0.0_WP,PREFRds)
    where (new_snow) REFR(:,:,month) = max(0.0_WP,PREFRns)
    where (mask(:,:,month))  REFR(:,:,month) = min(0.6_WP*snh_est/ml, REFR(:,:,month))
    if (debug_switch) then
      write(*,*) "snh_est",snh_est(debug_lon, debug_lat)
      write(*,*) "SNOW(:,:,month)",SNOW(debug_lon, debug_lat,month)
      write(*,*) "PREFRws",PREFRws(debug_lon, debug_lat)
      write(*,*) "PREFRds",PREFRds(debug_lon, debug_lat)
      write(*,*) "PREFRns",PREFRns(debug_lon, debug_lat)
      write(*,*) "REFR",REFR(debug_lon, debug_lat,month)
    end if

    ! SNH
    SNH(:,:,month)   = max(0.0_WP,(tmpSNH + (SNOW(:,:,month) + REFR(:,:,month) - MELT(:,:,month))*ml))
    if (debug_switch) then
      write(*,*) "SNH...."
      write(*,*) "SNH",SNH(debug_lon, debug_lat,month)
    end if

    ! Albedo ~ snow type
    where (wet_snow) A(:,:,month) = Aws
    where (dry_snow) A(:,:,month) = Adsc
    where (new_snow) A(:,:,month) = Ansc
    if (debug_switch) then
      write(*,*) "A",A(debug_lon, debug_lat,month)
    end if

    ! SMB
   where (mask(:,:,month))  SMB(:,:,month) = SNOW(:,:,month) - MELT(:,:,month) + REFR(:,:,month)
   if (debug_switch) then
     write(*,*) "SMB",SMB(debug_lon, debug_lat,month)
   end if

   if (month==9) then
     tmpSNH = max(0.0_WP,(SNH(:,:,month) - lastyear))
     SNH(:,:,month) = tmpSNH
     lastyear = tmpSNH
   else
     tmpSNH = SNH(:,:,month)
   end if

end do
END SUBROUTINE dEBM_core

END MODULE MOD_MAIN
