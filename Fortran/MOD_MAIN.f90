MODULE MOD_MAIN
  USE MOD_PRE
  USE MOD_DATA

  IMPLICIT NONE

contains

SUBROUTINE dEBM_sunny_hours_c
  !******************************************************************************
  ! calculates the time, the sun is above a certain elevation angle
  ! elevation needs to be a scalar
  ! lat may be scalar, vector or matrix
  ! declination needs to be scalar, or same size as lat

  real :: elev1
  real, allocatable, dimension(:,:,:,:) :: decl
  real, allocatable, dimension(:,:,:,:) :: sinphisind,cosphicosd
  real, allocatable, dimension(:,:,:,:) :: ha_0,ha_elev
  real, allocatable, dimension(:,:,:,:) :: sol_flux_fact_0,sol_flux_fact_e

  allocate (decl(xlen,ylen,mlen,nlen))
  allocate (sinphisind(xlen,ylen,mlen,nlen))
  allocate (cosphicosd(xlen,ylen,mlen,nlen))
  allocate (ha_0(xlen,ylen,mlen,nlen))
  allocate (ha_elev(xlen,ylen,mlen,nlen))
  allocate (sol_flux_fact_0(xlen,ylen,mlen,nlen))
  allocate (sol_flux_fact_e(xlen,ylen,mlen,nlen))

  decl = spread(spread(spread( obl*sin(pi/180*360/365*(days-79)), 1, xlen ), 2 ,ylen), 4 ,nlen)

  ! elev
  elev1            = elev*pi/180
  !
  ! !products of sin(lat)*sin(decl) and cos(lat)*cos(decl)  is needed more than once
  sinphisind       = sin(pi/180*lat)*sin(pi/180*decl)
  cosphicosd       = cos(pi/180*lat)*cos(pi/180*decl)

  ! ! hourangle of sun rise and factor which relates solar flux density and
  ! ! daily insolation SW
   ha_0             = real(acos(-sinphisind/cosphicosd))
   sol_flux_fact_0  = (sinphisind*ha_0+cosphicosd*sin(ha_0))/pi

  ! hourangle of sun rising above elev
  ha_elev           = real(acos((sin(elev1)-sinphisind)/cosphicosd))

  ! duration of melt period in hours
  hours            = ha_elev*24/pi

  ! proportion of daily insolation which is effective during melt period
  sol_flux_fact_e  = (sinphisind*ha_elev+cosphicosd*sin(ha_elev))/pi
  q                = sol_flux_fact_e/sol_flux_fact_0

END SUBROUTINE dEBM_sunny_hours_c





SUBROUTINE PDD4

  PDD=stddev/sqrt(2*pi)*exp(-(surface_temperature1**2)/(2*(stddev**2)))+surface_temperature1/2.*erfc(-surface_temperature1/(sqrt(2.)*stddev))

END SUBROUTINE PDD4




SUBROUTINE dEBM_main

  real, allocatable, dimension(:,:,:,:) :: shortwave_radiation_downward0,melt

  allocate (melt(xlen,ylen,mlen,nlen))
  allocate (shortwave_radiation_downward0(xlen,ylen,mlen,nlen))

  melt(:,:,:,:) = 0
  ! melt condition: ice is warm enough to warm to melting point during daytime
  where (surface_temperature1>Tmin) melt(:,:,:,:) = 1.

  ! ! calculate monthly melt rates
  ! ! dEBMmodel
  shortwave_radiation_downward0 = (1-albedo1)*shortwave_radiation_downward1*24*q
  where (shortwave_radiation_downward0 < 0) shortwave_radiation_downward0 = 0.
  dEBMmelt=(melt*shortwave_radiation_downward0+(c2+c1*PDD)*hours)/Lf*60*60

END SUBROUTINE dEBM_main

END MODULE MOD_MAIN
