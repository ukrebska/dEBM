MODULE MOD_OUTPUT
 USE MOD_PRE
 USE MOD_DATA
 USE MOD_MAIN

 IMPLICIT NONE

contains

SUBROUTINE OUTPUT_NETCDF
  USE MOD_PRE
  USE MOD_DATA
  USE MOD_MAIN
  IMPLICIT NONE

  include 'netcdf.inc'

  character (len = *), parameter :: filename_out = "surface_melt.nc"
  integer :: ncid, status

  integer, parameter :: NDIMS = 3
  integer :: NLONS, NLATS, NTIMS
  character (len = *), parameter :: LON_NAME = "x"
  character (len = *), parameter :: LAT_NAME = "y"
  character (len = *), parameter :: TIM_NAME = "t"
  integer :: lon_dimid,lat_dimid,  tim_dimid

  real, dimension(:), allocatable :: lons, lats, time
  integer :: lon_varid, lat_varid,  tim_varid

  character (len = *), parameter :: MELT_NAME="smelt"
  integer :: melt_varid

  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: MELT_UNITS = "mm/day"
  character (len = *), parameter :: LON_UNITS = "degrees_east"
  character (len = *), parameter :: LAT_UNITS = "degrees_north"
  character (len = *), parameter :: TIM_UNITS = "day since 1-1-1"

  real, dimension(:,:,:), allocatable :: melt_out

  NLONS = xlen
  NLATS = ylen
  NTIMS = tlen

  ! Allocate memory.
  allocate(melt_out(NLONS, NLATS, NTIMS))
  allocate(lons(NLONS))
  allocate(lats(NLATS))
  allocate(time(NTIMS))

  lons=x
  lats=y
  time=t
  melt_out= reshape( dEBMmelt, (/NLONS,NLATS,NTIMS/) )

  ! Create the file.
  status=nf_create(filename_out, NF_NETCDF4, ncid)
  if (status .ne. nf_noerr) call handle_err(status)
  print *,"creating ",trim(filename_out)," ..."
  !
  ! ! Define the dimensions.
  status=nf_def_dim(ncid, LAT_NAME, NLATS, lat_dimid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_def_dim(ncid, LON_NAME, NLONS, lon_dimid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_def_dim(ncid, TIM_NAME, NTIMS, tim_dimid)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_def_var(ncid, LAT_NAME, nf_real, 1, lat_dimid, lat_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_def_var(ncid, LON_NAME, nf_real, 1, lon_dimid, lon_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_def_var(ncid, TIM_NAME, nf_real, 1, tim_dimid, tim_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_put_att_text(ncid, lat_varid, UNITS, len(LAT_UNITS), LAT_UNITS)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_put_att_text(ncid, lon_varid, UNITS, len(LON_UNITS), LON_UNITS)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_put_att_text(ncid, tim_varid, UNITS, len(TIM_UNITS), TIM_UNITS)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_def_var(ncid, MELT_NAME, NF_real, NDIMS,  (/ lon_dimid, lat_dimid, tim_dimid /), melt_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, melt_varid, UNITS, len(MELT_UNITS), MELT_UNITS)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_enddef(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_put_var_real(ncid, lat_varid, lats)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_var_real(ncid, lon_varid, lons)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_var_real(ncid, tim_varid, time)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_put_var_real(ncid, melt_varid, melt_out)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  print *," **** Done. ****"
END SUBROUTINE OUTPUT_NETCDF

SUBROUTINE handle_err(errcode)
  implicit none
  include 'netcdf.inc'
  integer errcode

  print *, 'Error: ', nf_strerror(errcode)
  stop "stopped"
END SUBROUTINE handle_err

END MODULE MOD_OUTPUT
