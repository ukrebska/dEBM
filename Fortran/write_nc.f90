program OUTPUT_NETCDF
  IMPLICIT NONE

  include 'netcdf.inc'

  character (len = *), parameter :: FILE_OUT = "surface_melt.nc"
  integer :: ncid, status

  integer, parameter :: NDIMS = 3
  integer, parameter :: NLATS = 96, NLONS = 192, NTIMS = 12
  character (len = *), parameter :: LON_NAME = "x"
  character (len = *), parameter :: LAT_NAME = "y"
  character (len = *), parameter :: TIM_NAME = "t"
  integer :: lon_dimid,lat_dimid,  tim_dimid

  real :: lons(NLONS), lats(NLATS),  time(NTIMS)
  integer :: lon_varid, lat_varid,  tim_varid

  character (len = *), parameter :: MELT_NAME="smelt"
  integer :: melt_varid
  integer :: dimids(NDIMS)

  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: MELT_UNITS = "mm/day"
  character (len = *), parameter :: LON_UNITS = "degrees_east"
  character (len = *), parameter :: LAT_UNITS = "degrees_north"
  character (len = *), parameter :: TIM_UNITS = "day since 1-1-1"

  real, dimension(:,:,:), allocatable :: melt_out

  ! Allocate memory.
  allocate(melt_out(NLONS, NLATS, NTIMS))

  lons(NLONS)=1
  lats(NLATS)=1
  time(NTIMS)=1
  melt_out(NLONS, NLATS, NTIMS)=9

  ! Create the file.
  status=nf_create(FILE_OUT, nf_netcdf4, ncid)
  !
  ! ! Define the dimensions.
  status=nf_def_dim(ncid, LAT_NAME, NLATS, lat_dimid)
  status=nf_def_dim(ncid, LON_NAME, NLONS, lon_dimid)
  status=nf_def_dim(ncid, TIM_NAME, NTIMS, tim_dimid)

  status=nf_def_var(ncid, LAT_NAME, NF_DOUBLE, 1, lat_dimid, lat_varid)
  status=nf_def_var(ncid, LON_NAME, NF_DOUBLE, 1, lon_dimid, lon_varid)
  status=nf_def_var(ncid, TIM_NAME, NF_DOUBLE, 1, tim_dimid, tim_varid)

  status=nf_put_att_text(ncid, lat_varid, UNITS, len(LAT_UNITS), LAT_UNITS)
  status=nf_put_att_text(ncid, lon_varid, UNITS, len(LON_UNITS), LON_UNITS)
  status=nf_put_att_text(ncid, tim_varid, UNITS, len(TIM_UNITS), TIM_UNITS)

  status = nf_def_var(ncid, MELT_NAME, NF_DOUBLE, NDIMS,  (/ lon_dimid, lat_dimid, tim_dimid /), melt_varid)
  status = nf_put_att_text(ncid, melt_varid, UNITS, len(MELT_UNITS), MELT_UNITS)

  status=nf_enddef(ncid)

  status = nf_put_var_double(ncid, lat_varid, lats)
  status = nf_put_var_double(ncid, lon_varid, lons)
  status = nf_put_var_double(ncid, tim_varid, time)

  status = nf_put_var_double(ncid, melt_varid, melt_out)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)

  print *,"*** SUCCESS writing example file surface_melt.nc!"

end program OUTPUT_NETCDF


SUBROUTINE handle_err(errcode)
implicit none
include 'netcdf.inc'
integer errcode

print *, 'Error: ', nf_strerror(errcode)
stop 2
END SUBROUTINE handle_err
