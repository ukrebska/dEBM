PROGRAM dEBMmain
  !******************************************************************************
  ! please cite:
  ! Krebs-Kanzow, U., Gierz, P., and Lohmann, G.: Brief communication:
  ! An ice surface melt scheme including the diurnal cycle of solar radiation,
  ! The Cryosphere, 12, 3923–3930, https://doi.org/10.5194/tc-12-3923-2018, 2018.
  !******************************************************************************
  ! Main calculations for dEBM
  !  | MOD_PRE: reads initial prescribed parameters from namelists
  !  | MOD_DATA: reads the input atmospheric forcing
  !  | MOD_MAIN: calcualtes SMB
  !  | MOD_OUTPUT: writes output files in NC format
  !  | MOD_RESTART: writes restart files in NC format
  !
  ! Coded by Shan Xu
  !******************************************************************************

  USE MOD_PRE
  USE MOD_DATA
  USE MOD_MAIN
  USE MOD_OUTPUT
  USE MOD_RESTART

  implicit none

  integer :: n
  integer :: mth_str

  real(kind=WP), dimension(:,:), allocatable     :: snh_tmp, snh_lastyear
  real(kind=WP), dimension(:,:,:,:), allocatable ::  snow_height, surface_mass_balance, &
                                              &melt_rate, accmulation_rate, &
                                              &refreeze_rate, albedo,&
                                              &snow_amount, rain_rate
  real(kind=WP), dimension(:), allocatable :: summer_solar_density

  write(*,*) "================== Begin =================="

  write(*,*) "Starting calculating dEBM..."

  write(*,*) "Reading in parameters from namelist..."
  CALL read_namelist

  write(*,*) "Reading in atmosphere fields...."
  CALL get_init

  ! init snow height
  allocate (snh_tmp( xlen, ylen))
  allocate (snh_lastyear( xlen, ylen))

  ! init main SMB components
  allocate (snow_height( xlen, ylen, mlen, nlen))
  allocate (surface_mass_balance( xlen, ylen, mlen, nlen))
  allocate (melt_rate( xlen, ylen, mlen, nlen))
  allocate (accmulation_rate( xlen, ylen, mlen, nlen))
  allocate (refreeze_rate( xlen, ylen, mlen, nlen))
  allocate (albedo( xlen, ylen, mlen, nlen))
  allocate (snow_amount( xlen, ylen, mlen, nlen))
  allocate (rain_rate( xlen, ylen, mlen, nlen))
  allocate (summer_solar_density(nlen))
  snow_height(:,:,:,:)          = 0.0_WP
  surface_mass_balance(:,:,:,:) = 0.0_WP
  melt_rate(:,:,:,:)            = 0.0_WP
  accmulation_rate(:,:,:,:)     = 0.0_WP
  refreeze_rate(:,:,:,:)        = 0.0_WP
  albedo(:,:,:,:)               = 0.0_WP
  snow_amount(:,:,:,:)          = 0.0_WP
  rain_rate(:,:,:,:)            = 0.0_WP
  summer_solar_density(:)       = 0.0_WP

  ! Init snow height from restart or spin-up
  n=1
  if ((n==1) .AND. (lresume==.true.))then
    ! For a warm start, read the snow height from restart file
    write(*,*) "We are doing a restart."
    CALL read_restart(snh_tmp, snh_lastyear)
  else
    ! For a cold start, we'll do spin–up based on the first 15 months
    write(*,*) "We are doing a cold start."
    write(*,*) "The first 15 months is regarded as spin-up as default."
    write(*,*) "If you want to restart, plz prescribe 'lresume=.true.' in the namelist"
    ! First, we integrate from <Sep to Dec>; then use the forcing of the first year for another 12-month cycles.
    ! init snow height with no snow cover
    snh_tmp(:,:)      = 0.0_WP
    snh_lastyear(:,:) = 0.0_WP
    ! First, we calculate the first year from Oct to Dec
    mth_str = 9     ! start from Sep
    CALL dEBM_core(surface_temperature1(:,:,:,n), &
                    &shortwave_radiation_downward1(:,:,:,n), shortwave_radiation_TOA1(:,:,:,n), &
                    &emissivity1(:,:,:,n), cloud_cover1(:,:,:,n), precipitation1(:,:,:,n), snh_tmp, snh_lastyear, &
                    &lat(:,:,:,n), mask1, obliquity, mth_str, &
                    &snow_height(:,:,:,n), surface_mass_balance(:,:,:,n), melt_rate(:,:,:,n), refreeze_rate(:,:,:,n), albedo(:,:,:,n), &
                    &snow_amount(:,:,:,n), rain_rate(:,:,:,n), summer_solar_density(n))
    snh_tmp      = snow_height(:,:,12,1)
    snh_lastyear = snow_height(:,:,9,1)
    ! Second, we recalculate the first year from Jan to Dec
    mth_str = 1     ! start from Jan
    CALL dEBM_core(surface_temperature1(:,:,:,n), &
                    &shortwave_radiation_downward1(:,:,:,n), shortwave_radiation_TOA1(:,:,:,n), &
                    &emissivity1(:,:,:,n), cloud_cover1(:,:,:,n), precipitation1(:,:,:,n), snh_tmp, snh_lastyear, &
                    &lat(:,:,:,n), mask1, obliquity, mth_str, &
                    &snow_height(:,:,:,n), surface_mass_balance(:,:,:,n), melt_rate(:,:,:,n), refreeze_rate(:,:,:,n), albedo(:,:,:,n), &
                    &snow_amount(:,:,:,n), rain_rate(:,:,:,n), summer_solar_density(n))
    snh_tmp      = snow_height(:,:,12,1)
    snh_lastyear = snow_height(:,:,9,1)
  end if

  ! Actual simulation
  do n=1, nlen
    write(*,*) "calculates year",n
    mth_str = 1
    CALL dEBM_core(surface_temperature1(:,:,:,n), &
                    &shortwave_radiation_downward1(:,:,:,n), shortwave_radiation_TOA1(:,:,:,n), &
                    &emissivity1(:,:,:,n), cloud_cover1(:,:,:,n), precipitation1(:,:,:,n), snh_tmp, snh_lastyear, &
                    &lat(:,:,:,n), mask1, obliquity, mth_str, &
                    &snow_height(:,:,:,n), surface_mass_balance(:,:,:,n), melt_rate(:,:,:,n), refreeze_rate(:,:,:,n), albedo(:,:,:,n), &
                    &snow_amount(:,:,:,n), rain_rate(:,:,:,n), summer_solar_density(n))
    snh_tmp      = snow_height(:,:,12,n)
    snh_lastyear = snow_height(:,:,9,n)
  end do

  ! Debug
  if (debug_switch) then
    write(*,*) "SNH",snow_height(debug_lon, debug_lat, debug_mon, debug_year)
    write(*,*) "surface_mass_balance",surface_mass_balance(debug_lon, debug_lat, debug_mon, debug_year)
    write(*,*) "melt_rate",melt_rate(debug_lon, debug_lat, debug_mon, debug_year)
    write(*,*) "albedo",albedo(debug_lon, debug_lat, debug_mon, debug_year)
    write(*,*) "summer_solar_density",summer_solar_density(7)
  end if

  ! Write restart
  CALL write_restart(snh_tmp, snh_lastyear)

  ! Convert output data to expected units
  snow_height           = snow_height/1000.                         ! SNH: convert units from "mm" to "m"
  surface_mass_balance  = surface_mass_balance/(24.*60.*60.)        ! SMB: convert units from "mm day-1" to "kg m-2 second-1"
  melt_rate             = melt_rate/(24.*60.*60.)                   ! ME: convert units from "mm day-1" to "kg m-2 second-1"
  refreeze_rate         = refreeze_rate/(24.*60.*60.)               ! RZ: convert units from "mm day-1" to "kg m-2 second-1"
  snow_amount           = snow_amount/(24.*60.*60.)                 ! SNOW: convert units from "mm day-1" to "kg m-2 second-1"
  rain_rate             = rain_rate/(24.*60.*60.)                   ! RAIN: convert units from "mm day-1" to "kg m-2 second-1"
  ! Write_output
  CALL write_output(lon0, lat0, snow_height, surface_mass_balance, melt_rate,&
                      &refreeze_rate, albedo,&
                      &snow_amount, rain_rate)

  write(*,*) "================== Done =================="

END PROGRAM dEBMmain
