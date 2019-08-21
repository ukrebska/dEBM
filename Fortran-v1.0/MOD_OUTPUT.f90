MODULE MOD_OUTPUT
! ************************************************************************
! * This MODULE writes output files for the model                        *
! ************************************************************************
! * Variables includes: SNH, SMB, MELT, ACC, REFR, A                     *
! *   SNH           : snow height                                        *
! *                   units meter per month                              *
! *   SMB           : surface mass balance                               *
! *                   units kg m-2 second-1                              *
! *   MELT          : surface melt                                       *
! *                   units kg m-2 second-1                              *
! *   ACC           : accumulation                                       *
! *                   units kg m-2 second-1                              *
! *   REFR          : refreeze                                           *
! *                   units kg m-2 second-1                              *
! *   A             : albedo                                             *
! ************************************************************************

 USE MOD_PRE
 USE MOD_DATA
 USE MOD_MAIN

 IMPLICIT NONE

contains

  SUBROUTINE OUTPUT_NETCDF(lon_output, lat_output, snh_output, smb_output, melt_output, acc_output, refr_output, albedo_output)
    IMPLICIT NONE

    include 'netcdf.inc'

    real, intent(in), dimension(:,:) :: lon_output, lat_output
    real, intent(in), dimension(:,:,:,:) :: snh_output, smb_output, melt_output, acc_output, refr_output, albedo_output

    character (len = *), parameter :: filename_out = "surface_mass_balance.nc"
    integer :: ncid, status

    integer, parameter :: NDIMS = 3
    integer :: NLONS, NLATS, NTIMS
    integer :: m, n, start
    character (len = *), parameter :: LON_NAME = "x"
    character (len = *), parameter :: LAT_NAME = "y"
    character (len = *), parameter :: TIM_NAME = "time"
    integer :: lon_dimid, lat_dimid,  tim_dimid
    integer :: t_dimid, y_dimid, x_dimid
    integer :: t_varid, y_varid, x_varid

    real, dimension(:), allocatable :: x1, y1, t1
    real, dimension(:,:), allocatable :: lons, lats
    ! real, dimension(:) :: time

    integer :: tim_varid, lat_varid, lon_varid
    integer :: melt_varid, snh_varid, smb_varid,&
                &acc_varid, refr_varid, albedo_varid

    character (len = *), parameter :: MELT_NAME = "MELT"
    character (len = *), parameter :: SNH_NAME = "SNH"
    character (len = *), parameter :: SMB_NAME = "SMB"
    character (len = *), parameter :: ACC_NAME = "ACC"
    character (len = *), parameter :: REFR_NAME = "REFR"
    character (len = *), parameter :: Albedo_NAME = "A"

    character (len = *), parameter :: UNITS = "units"
    character (len = *), parameter :: LON_UNITS = "degrees_east"
    character (len = *), parameter :: LAT_UNITS = "degrees_north"
    ! character (len = *), parameter :: TIM_UNITS = "day since 1-1-1"
    character (len = *), parameter :: MELT_UNITS = "kg m-2 second-1"
    character (len = *), parameter :: SNH_UNITS = "m"
    character (len = *), parameter :: SMB_UNITS = "kg m-2 second-1"
    character (len = *), parameter :: ACC_UNITS = "kg m-2 second-1"
    character (len = *), parameter :: REFR_UNITS = "kg m-2 second-1"
    character (len = *), parameter :: Albedo_UNITS = ""

    ! write(*,*) "1"

    NLONS = xlen
    NLATS = ylen
    NTIMS = tlen

    ! Allocate memory.
    allocate(lons(NLONS, NLATS))
    allocate(lats(NLONS, NLATS))
    ! write(*,*) "2"

    ! get coordinate information from input files
    allocate(t1(tlen))
    allocate(y1(ylen))
    allocate(x1(xlen))
    ! write(*,*) shape(t1)
    ! write(*,*) shape(y1)
    ! write(*,*) shape(x1)
    CALL READ_COORDINATE_INFO(t1, y1, x1)
    ! write(*,*) "time is",t1
    ! write(*,*) "x is",x1
    ! write(*,*) "y is",y1
    ! write(*,*) "3"

    ! Create the file.
    status = nf_create(filename_out, nf_clobber, ncid)
    if (status .ne. nf_noerr) call handle_err(status)
    ! print *,"creating ",trim(filename_out)," ..."

    ! Define the dimensions.
    status = nf_def_dim(ncid, 'time', NF_UNLIMITED, t_dimid)
    if (status .ne. nf_noerr) call handle_err(status)
    ! write(*,*) "4.0"
    status = nf_def_dim(ncid, 'x', NLONS, x_dimid)
    if (status .ne. nf_noerr) call handle_err(status)
    ! write(*,*) "4.1"
    status = nf_def_dim(ncid, 'y', NLATS, y_dimid)
    if (status .ne. nf_noerr) call handle_err(status)
    ! write(*,*) "4.2"

    ! Define the dimension variables
    status = nf_def_var(ncid, 'time', NF_FLOAT, 1, t_dimid, t_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    ! write(*,*) "4.3"
    status = nf_def_var(ncid, 'x', NF_FLOAT, 1, x_dimid, x_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    ! write(*,*) "4.4"
    status = nf_def_var(ncid, 'y', NF_FLOAT, 1, y_dimid, y_varid)
    if (status .ne. nf_noerr) call handle_err(status)

    ! write(*,*) "4"

    ! ! Define the coordinate variables
    ! status = nf_def_var(ncid, 'lon', NF_FLOAT, 2, (/x_dimid, y_dimid/), lon_varid)
    ! if (status .ne. nf_noerr) call handle_err(status)
    ! status = nf_def_var(ncid, 'lat', NF_FLOAT, 2, (/x_dimid, y_dimid/), lat_varid)
    ! if (status .ne. nf_noerr) call handle_err(status)
    ! ! status=nf_def_var(ncid, TIM_NAME, nf_init, 1, tim_dimid, tim_varid)
    ! ! if (status .ne. nf_noerr) call handle_err(status)
    !
    ! write(*,*) "5"
    !
    ! ! Assign units attributes to coordinate variables.
    ! status = nf_put_att_text(ncid, lon_varid, UNITS, len(LON_UNITS), LON_UNITS)
    ! if (status .ne. nf_noerr) call handle_err(status)
    ! status = nf_put_att_text(ncid, lat_varid, UNITS, len(LAT_UNITS), LAT_UNITS)
    ! if (status .ne. nf_noerr) call handle_err(status)
    ! ! status=nf_put_att_text(ncid, tim_varid, UNITS, len(TIM_UNITS), TIM_UNITS)
    ! ! if (status .ne. nf_noerr) call handle_err(status)

    ! write(*,*) "6"

    ! snh
    status = nf_def_var(ncid, SNH_NAME, NF_FLOAT, NDIMS, (/x_dimid, y_dimid, t_dimid/), snh_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_put_att_text(ncid, snh_varid, UNITS, len(SNH_UNITS), SNH_UNITS)
    if (status .ne. nf_noerr) call handle_err(status)

    ! melt
    status = nf_def_var(ncid, MELT_NAME, NF_FLOAT, NDIMS, (/x_dimid, y_dimid, t_dimid/), melt_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_put_att_text(ncid, melt_varid, UNITS, len(MELT_UNITS), MELT_UNITS)
    if (status .ne. nf_noerr) call handle_err(status)

    ! smb
    status = nf_def_var(ncid, SMB_NAME, NF_FLOAT, NDIMS, (/x_dimid, y_dimid, t_dimid/), smb_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_put_att_text(ncid, smb_varid, UNITS, len(SMB_UNITS), SMB_UNITS)
    if (status .ne. nf_noerr) call handle_err(status)

    ! acc
    status = nf_def_var(ncid, ACC_NAME, NF_FLOAT, NDIMS, (/x_dimid, y_dimid, t_dimid/), acc_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_put_att_text(ncid, acc_varid, UNITS, len(ACC_UNITS), ACC_UNITS)
    if (status .ne. nf_noerr) call handle_err(status)

    ! write(*,*) "7"

    ! refr
    status = nf_def_var(ncid, REFR_NAME, NF_FLOAT, NDIMS, (/x_dimid, y_dimid, t_dimid/), refr_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_put_att_text(ncid, refr_varid, UNITS, len(REFR_UNITS), REFR_UNITS)
    if (status .ne. nf_noerr) call handle_err(status)

    ! albedo
    status = nf_def_var(ncid, Albedo_NAME, NF_FLOAT, NDIMS, (/x_dimid, y_dimid, t_dimid/), Albedo_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_put_att_text(ncid, Albedo_varid, UNITS, len(Albedo_UNITS), Albedo_UNITS)
    if (status .ne. nf_noerr) call handle_err(status)

    ! write(*,*) "8"

    ! Close define mode
    status = nf_enddef(ncid)
    if (status .ne. nf_noerr) call handle_err(status)

    ! write(*,*) "9"

    ! Write the data.
    ! write dimension information
    ! write into the output file
    status=nf_put_vara_real(ncid, t_varid, 1, tlen, t1)
    if (status .ne. nf_noerr) call handle_err(status)
    ! write(*,*) "t1 is", t
  ! write(*,*) "9.1"
    status=nf_put_vara_real(ncid, y_varid, 1, ylen, y1)
    if (status .ne. nf_noerr) call handle_err(status)
    ! write(*,*) "9.2"
    status=nf_put_vara_real(ncid, x_varid, 1, xlen, x1)
    if (status .ne. nf_noerr) call handle_err(status)


    ! ! Write the coordinate variable data() latitudes and longitudes) into the netCDF file.
    ! status = nf_put_var_real(ncid, lat_varid, lat_output)
    ! if (status .ne. nf_noerr) call handle_err(status)
    ! status = nf_put_var_real(ncid, lon_varid, lon_output)
    ! if (status .ne. nf_noerr) call handle_err(status)

    ! write(*,*) "10"

    ! write 3D variables
    start = 1
    do n = 1, nlen
      do m = 1, mlen

         ! snh
         status = nf_put_vara_real(ncid, snh_varid, (/1, 1, start/), (/NLONS, NLATS, 1/), snh_output(:, :, m, n))
         if (status .ne. nf_noerr) call handle_err(status)
         ! melt
         status = nf_put_vara_real(ncid, melt_varid, (/1, 1, start/), (/NLONS, NLATS, 1/), melt_output(:, :, m, n))
         if (status .ne. nf_noerr) call handle_err(status)
         ! smb
         status = nf_put_vara_real(ncid, smb_varid, (/1, 1, start/), (/NLONS, NLATS, 1/), smb_output(:, :, m, n))
         if (status .ne. nf_noerr) call handle_err(status)
         ! acc
         status = nf_put_vara_real(ncid, acc_varid, (/1, 1, start/), (/NLONS, NLATS, 1/), acc_output(:, :, m, n))
         if (status .ne. nf_noerr) call handle_err(status)
         ! refr
         status = nf_put_vara_real(ncid, refr_varid, (/1, 1, start/), (/NLONS, NLATS, 1/), refr_output(:, :, m, n))
         if (status .ne. nf_noerr) call handle_err(status)
         ! albedo
         status = nf_put_vara_real(ncid, albedo_varid, (/1, 1, start/), (/NLONS, NLATS, 1/), albedo_output(:, :, m, n))
         if (status .ne. nf_noerr) call handle_err(status)

         start = start + 1

       end do
    end do

    ! write(*,*) "11"

    ! close file
    status=nf_close(ncid)
    if (status .ne. nf_noerr) call handle_err(status)

  END SUBROUTINE OUTPUT_NETCDF

  SUBROUTINE READ_COORDINATE_INFO(t0, y0, x0)
    implicit NONE
    include 'netcdf.inc'
    integer :: ncid0, status0

    character (len = *), parameter :: x_NAME = "x"
    character (len = *), parameter :: y_NAME = "y"
    character (len = *), parameter :: t_NAME = "time"
    integer :: tid, yid, xid
    real, intent(out), dimension(:)  :: t0, y0, x0

    ! open file and reading coordinating info from input files
    status0=nf_open(trim(filename_in), nf_nowrite, ncid0)
    if (status0 .ne. nf_noerr) call handle_err(status0)

    ! time
    status0=nf_inq_varid(ncid0, t_NAME, tid)
    status0=nf_get_vara_real(ncid0,tid,1,tlen,t0)
    ! write(*,*) "t is ",t0
    ! write(*,*) "ylen is ",ylen
    ! write(*,*) "LAT_NAME is ",LAT_NAME

    ! y
    status0=nf_inq_varid(ncid0, y_NAME, yid)
    ! write(*,*) "LAT_NAME is "
    !
    status0=nf_get_vara_real(ncid0,yid,1,ylen,y0)
    ! write(*,*) "y is "

    ! x
    status0=nf_inq_varid(ncid0, x_NAME, xid)
    status0=nf_get_vara_real(ncid0,xid,1,xlen,x0)
    ! write(*,*) "x is "

    ! close file
    status0=nf_close(ncid0)
    if (status0 .ne. nf_noerr) call handle_err(status0)

  END SUBROUTINE READ_COORDINATE_INFO


  SUBROUTINE handle_err(errcode)
    implicit none
    include 'netcdf.inc'
    integer errcode

    print *, 'Error: ', nf_strerror(errcode)
    stop "stopped"
  END SUBROUTINE handle_err

END MODULE MOD_OUTPUT
