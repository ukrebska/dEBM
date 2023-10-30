MODULE MOD_RESTART
  ! ************************************************************************
  ! * This MODULE write and read restart files                             *
  ! ************************************************************************
  ! * MOD_RESTART                                                          *
  ! *         | write_restart                                              *
  ! *         | read_restart                                               *
  ! ************************************************************************

 USE MOD_PRE
 USE MOD_DATA
 USE MOD_MAIN
 implicit none

 character(len = *), parameter :: restart_filename = "restart_debm.nc"
 character(len = *), parameter :: tmpSNH_NAME = "tmpSNH"
 character(len = *), parameter :: lastyear_NAME = "lastyear"

contains


SUBROUTINE write_restart(tmpSNH, lastyear)
  ! SUBROUTINE write_restart(tempm, swdm, swd_TOAm, emiss, clcov, ppm, tmpSNH, lastyear, latm, mask)
  ! TODO：currently we only keep the minmal variables
  ! ************************************************************************
  ! * write_restart produces restart files                                 *
  ! *                                                                      *
  ! * Variables includes: tmpSNH, lastyear                                 *
  ! *   tmpSNH : snow height of September                                  *
  ! *   lastyear : ssnow height of December                                *
  ! ************************************************************************

    implicit none
    include 'netcdf.inc'

    real(kind=WP), intent(in), dimension(:,:)   :: tmpSNH, lastyear ! snow height of December and September will be saved

    integer :: NLONS, NLATS, NTIMS
    integer :: m, n, start
    integer :: y_dimid, x_dimid
    integer :: y_varid, x_varid
    integer :: tmpSNH_varid, lastyear_varid

    character(len = *), parameter :: UNITS = "units"
    character(len = *), parameter :: tmpSNH_UNITS = ""
    character(len = *), parameter :: lastyear_UNITS = ""

    NLONS = xlen
    NLATS = ylen

    ! Create the file.
    status = nf_create(restart_filename, nf_clobber, ncid)
    if (status .ne. nf_noerr) call handle_err(status)
    if (debug_switch) then
      write(*,*) "write_restart: survive creating ",trim(restart_filename)," ..."
    end if

    ! Define the dimensions.
    status = nf_def_dim(ncid, 'x', NLONS, x_dimid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_def_dim(ncid, 'y', NLATS, y_dimid)
    if (status .ne. nf_noerr) call handle_err(status)
    if (debug_switch) then
      write(*,*) "write_restart: survive define dimension"
    end if

    ! Define the dimension variables
    status = nf_def_var(ncid, 'x', NF_FLOAT, 1, x_dimid, x_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_def_var(ncid, 'y', NF_FLOAT, 1, y_dimid, y_varid)
    if (status .ne. nf_noerr) call handle_err(status)

    ! Define the coordinate variables
    status = nf_def_var(ncid, tmpSNH_NAME, NF_FLOAT, 2, (/x_dimid, y_dimid/), tmpSNH_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_def_var(ncid, lastyear_NAME, NF_FLOAT, 2, (/x_dimid, y_dimid/), lastyear_varid)
    if (status .ne. nf_noerr) call handle_err(status)

    ! ! Assign units attributes to coordinate variables.
    ! status = nf_put_att_text(ncid, tmpSNH_varid, UNITS, len(tmpSNH_UNITS), tmpSNH_UNITS)
    ! if (status .ne. nf_noerr) call handle_err(status)
    ! status = nf_put_att_text(ncid, lastyear_varid, UNITS, len(lastyear_UNITS), lastyear_UNITS)
    ! if (status .ne. nf_noerr) call handle_err(status)

    ! Close define mode
    status = nf_enddef(ncid)
    if (status .ne. nf_noerr) call handle_err(status)

    ! Write the data.
    status=nf_put_vara_double(ncid, y_varid, 1, ylen, y1)
    if (status .ne. nf_noerr) call handle_err(status)
    status=nf_put_vara_double(ncid, x_varid, 1, xlen, x1)
    if (status .ne. nf_noerr) call handle_err(status)
    if (debug_switch) then
      write(*,*) "write_restart: survive write dimension information"
    end if

    ! Write the coordinate variable data tmpSNH and lastyear into the netCDF file.
    status = nf_put_vara_double(ncid, tmpSNH_varid, (/1, 1/), (/NLONS, NLATS/), tmpSNH)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_put_vara_double(ncid, lastyear_varid, (/1, 1/), (/NLONS, NLATS/), lastyear)
    if (status .ne. nf_noerr) call handle_err(status)
    if (debug_switch) then
      write(*,*) "write_restart: survive write all variables"
    end if

    ! close file
    status=nf_close(ncid)
    if (status .ne. nf_noerr) call handle_err(status)

  END SUBROUTINE write_restart


  SUBROUTINE read_restart(tmpSNHr, lastyearr)
    ! ************************************************************************
    ! * read_restart produces restart files                                  *
    ! *                                                                      *
    ! * Variables includes: tmpSNH, lastyear                                 *
    ! *   tmpSNH : snow height of September                                  *
    ! *   lastyear : ssnow height of December                                *
    ! ************************************************************************

      implicit none
      include 'netcdf.inc'

      integer :: i,j
      integer :: y_dimid, x_dimid
      integer :: y_varid, x_varid
      integer :: xlen0, ylen0
      real(kind=WP), dimension(:), allocatable :: x0, y0
      real(kind=WP), intent(out), dimension(:,:), allocatable :: tmpSNHr, lastyearr ! snow height of December and September will be saved

      ! open file.
      status = nf_open(trim(restart_filename), nf_nowrite, ncid)
      if (status .ne. nf_noerr) then
        write(*,*) "cannot oepn ",trim(restart_filename)," ..."
        call handle_err(status)
      end if

      ! get the dimensions to check if restart files matches input
      ! x
      status = nf_inq_dimid(ncid, 'x', x_dimid)
      status = nf_inq_dimlen(ncid, x_dimid, xlen0)
      status = nf_inq_varid(ncid, 'x', x_varid)
      allocate(x0(xlen0))
      status = nf_get_vara_double(ncid, x_varid, 1, xlen0, x0)
      if (status .ne. nf_noerr) call handle_err(status)
      ! y
      status = nf_inq_dimid(ncid, 'y', y_dimid)
      status = nf_inq_dimlen(ncid, y_dimid, ylen0)
      status = nf_inq_varid(ncid, 'y', y_varid)
      allocate(y0(ylen0))
      status = nf_get_vara_double(ncid, y_varid, 1, ylen0, y0)
      if (status .ne. nf_noerr) call handle_err(status)
      ! check if restart files matches input
      if ((xlen0 /= xlen) .OR. (ylen0 /= ylen)) then
        write(*,*) "restart files didn't match input files"
        write(*,*) "mesh grid mismatch"
        write(*,*) " input mesh is ", xlen, "x", ylen
        write(*,*) " restart mesh  is ", xlen, "x", ylen
        stop "error detected!"
      else
        do i = 1, xlen0
          do j = 1, ylen0
            if ((x0(i) /= x1(i)) .OR. (y0(j) /= y1(j))) then
              write(*,*) "restart files didn't match input files"
              write(*,*) "mesh grid mismatch", i, j
              write(*,*) " input lon is ", x0(i), "lat is", y0(j)
              write(*,*) " input lon is ", x1(i), "lat is", y1(j)
              stop "error detected!"
            end if
          end do
        end do
      end if

      ! get tmpSNH
      allocate (tmpSNHr(xlen0,ylen0))
      status = nf_inq_varid(ncid, tmpSNH_NAME, varid)
      if (status .ne. nf_noerr) call handle_err(status)
      status = nf_get_vara_double(ncid, varid, (/1, 1/), (/xlen0, ylen0/), tmpSNHr)
      if (status .ne. nf_noerr) call handle_err(status)
      if (debug_switch) then
        write(*,*) "read_restart: survive reading ",tmpSNH_NAME
        write(*,*) "tmpSNH", tmpSNHr(debug_lon, debug_lat)
      end if

      ! get lastyear
      allocate (lastyearr(xlen0,ylen0))
      status = nf_inq_varid(ncid, lastyear_NAME, varid)
      status = nf_get_vara_double(ncid, varid, (/1, 1/), (/xlen0, ylen0/), lastyearr)
      if (status .ne. nf_noerr) call handle_err(status)
      if (debug_switch) then
        write(*,*) "read_restart: survive reading ",lastyear_NAME
        write(*,*) "lastyear", lastyearr(debug_lon, debug_lat)
      end if

      ! close file
      status=nf_close(ncid)
      if (status .ne. nf_noerr) call handle_err(status)

    END SUBROUTINE read_restart

END MODULE MOD_RESTART
