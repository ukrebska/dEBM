MODULE MOD_DATA
  USE MOD_PRE
  IMPLICIT NONE
  ! input variables
  real, allocatable :: shortwave_radiation_downward(:,:,:),albedo(:,:,:),surface_temperature(:,:,:)
  real, allocatable :: shortwave_radiation_downward1(:,:,:,:),albedo1(:,:,:,:),surface_temperature1(:,:,:,:),PDD(:,:,:,:)
  real, allocatable :: lat(:,:,:,:),hours(:,:,:,:),q(:,:,:,:),dEBMmelt(:,:,:,:)
  real, allocatable :: x(:),y(:),t(:)
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

  SUBROUTINE VAR_ALLOCATE    !
    allocate (shortwave_radiation_downward(xlen,ylen,tlen))
    allocate (surface_temperature(xlen,ylen,tlen))
    allocate (albedo(xlen,ylen,tlen))
    allocate (x(xlen))
    allocate (y(ylen))
    allocate (t(tlen))
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

  SUBROUTINE READ_NETCDF

    integer :: status, ncid, varid

    include 'netcdf.inc'

    status=nf_open(trim(filename_in), nf_nowrite, ncid)
    !
    status=nf_inq_varid(ncid, shortwave_radiation_downward_varname, varid)
    status=nf_get_var_real(ncid, varid, shortwave_radiation_downward)
    shortwave_radiation_downward1 = reshape( shortwave_radiation_downward, (/xlen,ylen,mlen,nlen/) )

    status=nf_inq_varid(ncid, temperature_varname, varid)
    status=nf_get_var_real(ncid, varid, surface_temperature)
    surface_temperature=surface_temperature-273.15
    surface_temperature1 = reshape( surface_temperature, (/xlen,ylen,mlen,nlen/) )

    status=nf_inq_varid(ncid, albedo_varname, varid)
    status=nf_get_var_real(ncid, varid, albedo)
    albedo1 = reshape( albedo, (/xlen,ylen,mlen,nlen/) )

    status=nf_inq_varid(ncid, longitude_varname, varid)  !inquire var
    status=nf_get_var_real(ncid,varid,x)  ! to get x

    status=nf_inq_varid(ncid, latitude_varname, varid)  !inquire var
    status=nf_get_var_real(ncid,varid,y)
    lat = spread(spread(spread( y, 1, xlen ), 3 ,mlen), 4 ,nlen)

    status=nf_inq_varid(ncid, time_varname, varid)  !inquire var
    status=nf_get_var_real(ncid,varid,t)

    status=nf_close(ncid)

    END SUBROUTINE

END MODULE MOD_DATA
