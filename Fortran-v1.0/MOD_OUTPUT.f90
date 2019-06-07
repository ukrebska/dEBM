MODULE MOD_OUTPUT
! ************************************************************************
! * This MODULE writes output files for the model                        *
! ************************************************************************
! * Variables includes: SNH, SMB, MELT, ACC, REFR, A                     *
! *   SNH           : snow height                                        *
! *   SMB           : surface mass balance                               *
! *   MELT          : surface melt                                       *
! *   ACC           : accumulation                                       *
! *   REFR          : refreeze                                           *
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

    real, dimension(:,:), allocatable :: lons, lats
    ! real, dimension(:) :: time

    integer :: lon_varid, lat_varid
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
    character (len = *), parameter :: MELT_UNITS = "mm/month"
    character (len = *), parameter :: SNH_UNITS = "mm/month"
    character (len = *), parameter :: SMB_UNITS = "mm/month"
    character (len = *), parameter :: ACC_UNITS = "mm/month"
    character (len = *), parameter :: REFR_UNITS = "mm/month"
    character (len = *), parameter :: Albedo_UNITS = " "

    write(*,*) "1"

    NLONS = xlen
    NLATS = ylen
    NTIMS = tlen

    write(*,*) "2"

    ! Allocate memory.
    allocate(lons(NLONS, NLATS))
    allocate(lats(NLONS, NLATS))

    write(*,*) "3"

    ! Create the file.
    status = nf_create(filename_out, nf_clobber, ncid)
    if (status .ne. nf_noerr) call handle_err(status)
    print *,"creating ",trim(filename_out)," ..."

    ! Define the dimensions.
    status = nf_def_dim(ncid, TIM_NAME, NF_UNLIMITED, tim_dimid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_def_dim(ncid, LON_NAME, NLONS, lon_dimid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_def_dim(ncid, LAT_NAME, NLATS, lat_dimid)
    if (status .ne. nf_noerr) call handle_err(status)
    !status = nf_def_dim(ncid, TIM_NAME, NF_UNLIMITED, tim_dimid)
    !if (status .ne. nf_noerr) call handle_err(status)

    write(*,*) "4"

    ! Define the coordinate variables
    status = nf_def_var(ncid, 'lon', nf_real, 2, (/lon_dimid, lat_dimid/), lon_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_def_var(ncid, 'lat', nf_real, 2, (/lon_dimid, lat_dimid/), lat_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    ! status=nf_def_var(ncid, TIM_NAME, nf_init, 1, tim_dimid, tim_varid)
    ! if (status .ne. nf_noerr) call handle_err(status)

    write(*,*) "5"

    ! Assign units attributes to coordinate variables.
    status = nf_put_att_text(ncid, lon_varid, UNITS, len(LON_UNITS), LON_UNITS)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_put_att_text(ncid, lat_varid, UNITS, len(LAT_UNITS), LAT_UNITS)
    if (status .ne. nf_noerr) call handle_err(status)
    ! status=nf_put_att_text(ncid, tim_varid, UNITS, len(TIM_UNITS), TIM_UNITS)
    ! if (status .ne. nf_noerr) call handle_err(status)

    write(*,*) "6"

    ! snh
    status = nf_def_var(ncid, SNH_NAME, NF_real, NDIMS, (/lon_dimid, lat_dimid, tim_dimid/), snh_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_put_att_text(ncid, snh_varid, UNITS, len(SNH_UNITS), SNH_UNITS)
    if (status .ne. nf_noerr) call handle_err(status)

    ! melt
    status = nf_def_var(ncid, MELT_NAME, NF_real, NDIMS, (/lon_dimid, lat_dimid, tim_dimid/), melt_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_put_att_text(ncid, melt_varid, UNITS, len(MELT_UNITS), MELT_UNITS)
    if (status .ne. nf_noerr) call handle_err(status)

    ! smb
    status = nf_def_var(ncid, SMB_NAME, NF_real, NDIMS, (/lon_dimid, lat_dimid, tim_dimid/), smb_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_put_att_text(ncid, smb_varid, UNITS, len(SMB_UNITS), SMB_UNITS)
    if (status .ne. nf_noerr) call handle_err(status)

    ! acc
    status = nf_def_var(ncid, ACC_NAME, NF_real, NDIMS, (/lon_dimid, lat_dimid, tim_dimid/), acc_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_put_att_text(ncid, acc_varid, UNITS, len(ACC_UNITS), ACC_UNITS)
    if (status .ne. nf_noerr) call handle_err(status)

    write(*,*) "7"

    ! refr
    status = nf_def_var(ncid, REFR_NAME, NF_real, NDIMS, (/lon_dimid, lat_dimid, tim_dimid/), refr_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_put_att_text(ncid, refr_varid, UNITS, len(REFR_UNITS), REFR_UNITS)
    if (status .ne. nf_noerr) call handle_err(status)

    ! albedo
    status = nf_def_var(ncid, Albedo_NAME, NF_real, NDIMS, (/lon_dimid, lat_dimid, tim_dimid/), Albedo_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_put_att_text(ncid, Albedo_varid, UNITS, len(Albedo_UNITS), Albedo_UNITS)
    if (status .ne. nf_noerr) call handle_err(status)

    write(*,*) "8"

    ! Close define mode
    status = nf_enddef(ncid)
    if (status .ne. nf_noerr) call handle_err(status)

    write(*,*) "9"

    ! Write the coordinate variable data() latitudes and longitudes) into the netCDF file.
    status = nf_put_var_real(ncid, lat_varid, lat_output)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_put_var_real(ncid, lon_varid, lon_output)
    if (status .ne. nf_noerr) call handle_err(status)
    ! status = nf_put_var(ncid, tim_varid, t)
    ! if (status .ne. nf_noerr) call handle_err(status)

    write(*,*) "10"

    ! Write the data.
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

    write(*,*) "11"

    ! close file
    status=nf_close(ncid)
    if (status .ne. nf_noerr) call handle_err(status)

  END SUBROUTINE OUTPUT_NETCDF


  SUBROUTINE handle_err(errcode)
    implicit none
    include 'netcdf.inc'
    integer errcode

    print *, 'Error: ', nf_strerror(errcode)
    stop "stopped"
  END SUBROUTINE handle_err

END MODULE MOD_OUTPUT
