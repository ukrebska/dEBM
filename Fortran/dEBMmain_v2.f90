MODULE PRE
  integer       :: xlen, ylen,tlen
  character*300 :: filename
  character*20  :: albedo_varname,temperature_varname,&
                     &shortwave_radiation_downward_varname,longitude_varname,latitude_varname
  real          :: stddev
contains
  SUBROUTINE read_namelist
    namelist /na/ xlen, ylen,tlen,&
                    &filename,albedo_varname,temperature_varname,shortwave_radiation_downward_varname,&
                    &longitude_varname,latitude_varname,stddev
    open(10, file="./namelist_pre_new.txt", status='old' )
    read (10, nml=na)
    ! write(*,*) filename,albedo_varname,temperature_varname,shortwave_radiation_downward_varname
    close (10)
 END SUBROUTINE read_namelist
END MODULE PRE






MODULE DATA_MODULE
  USE PRE
  IMPLICIT NONE
  ! input variables
  real, allocatable :: shortwave_radiation_downward(:,:,:),albedo(:,:,:),surface_temperature(:,:,:)
  real, allocatable :: shortwave_radiation_downward1(:,:,:,:),albedo1(:,:,:,:),surface_temperature1(:,:,:,:),PDD(:,:,:,:)
  real, allocatable :: lat(:,:,:,:),hours(:,:,:,:),q(:,:,:,:),dEBMmelt(:,:,:,:)
  real, allocatable :: x(:),y(:)
  ! real, allocatable :: lat(:,:),lon(:,:)
  integer :: xdim,ydim,tdim
  ! character*300 :: file_in
  ! real, dimension(:,:,:,:), allocatable :: shortwave_radiation_downward,A,PDD,dEBMmelt,hours,q,melt
  ! real, dimension(:,:),  allocatable :: lat
  ! integer :: xdim, ydim, tdim, nyear
  ! parameters defined for dEBM
  real, parameter :: epsa   = .76     ! atm. emissivity
  real, parameter :: epsi   = .95     ! emissivity of ice
  real, parameter :: beta   = 10      ! turbulent heat transfer coeff
  real, parameter :: bolz   = 5.67e-8 ! Stefan-Boltzmann constant
  real, parameter :: T0     = 273.15  ! melt point in K
  real, parameter :: Tmin   = -6.5    ! background melt condition
  real, parameter :: elev   = 17.5
  real, parameter :: obl    = 23.44
  integer, dimension(12) :: mth = (/1,2,3,4,5,6,7,8,9,10,11,12/)
  integer, dimension(12) :: days = (/16,45,75,105,136,166,197,228,258,289,319,350/)
  real, parameter :: Lf     = 3.34e5         ! Lf is latent heat of fusion
  real, parameter :: pi     = 3.1415
  integer, parameter :: mlen   = 12
  integer :: nlen
  real :: c1,c2

contains
  SUBROUTINE CALC_PARAMETER(c1, c2)
    real, intent(out) :: c1, c2
    c1     = (epsa*4*bolz*(T0**3)+beta)
    c2     = (-epsi+epsa*epsi)*bolz*(T0**4)
  END SUBROUTINE CALC_PARAMETER

  SUBROUTINE VAR_ALLOCATE
    !
    allocate (shortwave_radiation_downward(xlen,ylen,tlen))
    allocate (surface_temperature(xlen,ylen,tlen))
    allocate (albedo(xlen,ylen,tlen))
    allocate (x(xlen))
    allocate (y(ylen))
    !
    nlen = tlen/12
    allocate (lat(xlen,ylen,mlen,nlen))
    allocate (shortwave_radiation_downward1(xlen,ylen,mlen,nlen))
    allocate (surface_temperature1(xlen,ylen,mlen,nlen))
    allocate (albedo1(xlen,ylen,mlen,nlen))
    allocate (PDD(xlen,ylen,mlen,nlen))
    allocate (hours(xlen,ylen,mlen,nlen))
    allocate (q(xlen,ylen,mlen,nlen))
    allocate (dEBMmelt(xlen,ylen,mlen,nlen))

  END SUBROUTINE VAR_ALLOCATE

END MODULE DATA_MODULE





PROGRAM TEST
  USE PRE
  USE DATA_MODULE

  IMPLICIT NONE

  include 'netcdf.inc'

  integer :: status,ncid,varid

  CALL read_namelist
  CALL CALC_PARAMETER(c1, c2)
  CALL VAR_ALLOCATE

  ! filename='/Users/xushan/Downloads/dEBM-master/pi477_u_MM_2601_270001.01_echam2.nc'
  status=nf_open(trim(filename), nf_nowrite, ncid)
  !
  status=nf_inq_varid(ncid, shortwave_radiation_downward_varname, varid)
  status=nf_get_var_real(ncid, varid, shortwave_radiation_downward)
  shortwave_radiation_downward1 = reshape( shortwave_radiation_downward, (/xlen,ylen,mlen,nlen/) )
  ! print *,shortwave_radiation_downward(159,45,1)

  status=nf_inq_varid(ncid, temperature_varname, varid)
  status=nf_get_var_real(ncid, varid, surface_temperature)
  surface_temperature=surface_temperature-273.15
  surface_temperature1 = reshape( surface_temperature, (/xlen,ylen,mlen,nlen/) )
  ! print *,shape(surface_temperature)
  ! print *,surface_temperature1(159,45,12,:)
  ! print *,surface_temperature1(159,45,:,2)

  status=nf_inq_varid(ncid, albedo_varname, varid)
  status=nf_get_var_real(ncid, varid, albedo)
  albedo1 = reshape( albedo, (/xlen,ylen,mlen,nlen/) )
  ! print *,A(159,45,1)

  status=nf_inq_varid(ncid, longitude_varname, varid)  !inquire var
  status=nf_get_var_real(ncid,varid,x)  ! to get x

  status=nf_inq_varid(ncid, latitude_varname, varid)  !inquire var
  status=nf_get_var_real(ncid,varid,y)
  lat = spread(spread(spread( y, 1, xlen ), 3 ,mlen), 4 ,nlen)
  ! print *,"lat dim is", shape(lat)
  ! print *,lat(3,:)
  ! print *,lat(:,2)

  status=nf_close(ncid)

  ! generate parameters that depend on latitude and month
  CALL dEBM_sunny_hours_c
  print *,"hours  is",hours(159,45,:,2)
  print *,"q  is",q(159,45,:,2)

  ! calculte PDD
  CALL PDD4
  ! print *,shape(PDD)

  CALL dEBM_main
  ! print *,shape(dEBMmelt)
  print *," dEBMmelt is ",dEBMmelt(159,45,:,2)

END PROGRAM TEST





SUBROUTINE dEBM_main
  USE PRE
  USE DATA_MODULE

  implicit none

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




SUBROUTINE dEBM_sunny_hours_c
  !******************************************************************************
  ! calculates the time, the sun is above a certain elevation angle
  ! elevation needs to be a scalar
  ! lat may be scalar, vector or matrix
  ! declination needs to be scalar, or same size as lat
  USE PRE
  USE DATA_MODULE

  implicit none

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
  ! print *,"decl dim is", shape(decl)

  ! elev
  elev1            = elev*pi/180
  !
  ! !products of sin(lat)*sin(decl) and cos(lat)*cos(decl)  is needed more than once
  sinphisind       = sin(pi/180*lat)*sin(pi/180*decl)
  cosphicosd       = cos(pi/180*lat)*cos(pi/180*decl)
  ! print *,"sinphisind",sinphisind(155,24,12,1)
  ! print *,cosphicosd(155,24,12,1)

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
 USE PRE
 USE DATA_MODULE

  PDD=stddev/sqrt(2*pi)*exp(-(surface_temperature1**2)/(2*(stddev**2)))+surface_temperature1/2.*erfc(-surface_temperature1/(sqrt(2.)*stddev))

END SUBROUTINE PDD4
