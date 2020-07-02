PROGRAM MAIN

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

  write(*,*) "********* Begin *********"

  write(*,*) "Starting calculating dEBM..."

  write(*,*) "Reading from namelist..."
  CALL read_namelist

  write(*,*) "Reading in atmosphere fields...."
  CALL get_init

  allocate (snh_tmp( xlen, ylen))
  allocate (snh_lastyear( xlen, ylen))
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

  ! ************************************************************************
  ! the first year starts from September
  n=1
  if ((n==1) .AND. (lresume==.true.))then
    write(*,*) "We are doing a restart"
    CALL read_restart(snh_tmp, snh_lastyear)
  else
    snh_tmp(:,:)      = 0.0_WP
    snh_lastyear(:,:) = 0.0_WP
  end if
  mth_str = 9
  CALL dEBM_core(surface_temperature1(:,:,:,n), &
                  &shortwave_radiation_downward1(:,:,:,n), shortwave_radiation_TOA1(:,:,:,n), &
                  &emissivity1(:,:,:,n), cloud_cover1(:,:,:,n), precipitation1(:,:,:,n), snh_tmp, snh_lastyear, &
                  &lat(:,:,:,n), mask1, obliquity, mth_str, &
                  &snow_height(:,:,:,n), surface_mass_balance(:,:,:,n), melt_rate(:,:,:,n), refreeze_rate(:,:,:,n), albedo(:,:,:,n), &
                  &snow_amount(:,:,:,n), rain_rate(:,:,:,n), summer_solar_density(n))
  snh_tmp      = snow_height(:,:,12,1)
  snh_lastyear = snow_height(:,:,9,1)
  ! ************************************************************************
  ! recalculate the first year
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
  ! debug
  if (debug_switch) then
    write(*,*) "SNH",snow_height(debug_lon, debug_lat, debug_mon, debug_year)
    write(*,*) "surface_mass_balance",surface_mass_balance(debug_lon, debug_lat, debug_mon, debug_year)
    write(*,*) "melt_rate",melt_rate(debug_lon, debug_lat, debug_mon, debug_year)
    write(*,*) "albedo",albedo(debug_lon, debug_lat, debug_mon, debug_year)
    write(*,*) "summer_solar_density",summer_solar_density(1)
  end if
  ! ************************************************************************
  ! write_restart
  CALL write_restart(snh_tmp, snh_lastyear)
  ! n=nlen
  ! CALL write_restart(surface_temperature1(:,:,12,n), &
  !                 &shortwave_radiation_downward1(:,:,12,n), shortwave_radiation_TOA1(:,:,12,n), &
  !                 &emissivity1(:,:,12,n), cloud_cover1(:,:,12,n), precipitation1(:,:,12,n), snh_tmp, snh_lastyear, &
  !                 &lat(:,:,12,n), mask1, &
  !                 &snow_height(:,:,12,n), surface_mass_balance(:,:,12,n), melt_rate(:,:,12,n), refreeze_rate(:,:,12,n), albedo(:,:,12,n), &
  !                 &snow_amount(:,:,12,n), rain_rate(:,:,12,n), summer_solar_density(n))
  ! ************************************************************************
  ! write_output
  ! convert output data to expected units
  snow_height  = snow_height/1000.                ! SNH: convert units from "mm" to "m"
  surface_mass_balance  = surface_mass_balance/(24.*60.*60.)        ! SMB: convert units from "mm day-1" to "kg m-2 second-1"
  melt_rate = melt_rate/(24.*60.*60.)       ! MELT: convert units from "mm day-1" to "kg m-2 second-1"
  refreeze_rate = refreeze_rate/(24.*60.*60.)       ! REFR: convert units from "mm day-1" to "kg m-2 second-1"
  snow_amount = snow_amount/(24.*60.*60.)       ! REFR: convert units from "mm day-1" to "kg m-2 second-1"
  rain_rate = rain_rate/(24.*60.*60.)       ! REFR: convert units from "mm day-1" to "kg m-2 second-1"
  CALL write_output(lon0, lat0, snow_height, surface_mass_balance, melt_rate,&
                      &refreeze_rate, albedo,&
                      &snow_amount, rain_rate)

  write(*,*) "********* Done. *********"

END PROGRAM MAIN
