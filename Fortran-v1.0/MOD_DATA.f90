MODULE MOD_DATA
! ************************************************************************
! * This MODULE prepares for parameters and input variables              *
! ************************************************************************

  USE MOD_PRE
  IMPLICIT NONE

  ! Define variables
  real, dimension(:,:,:), allocatable   :: shortwave_radiation_downward,&
                                            &precipitation,&
                                            &surface_temperature,&
                                            &mask
  real, dimension(:,:,:,:), allocatable :: lat, shortwave_radiation_downward1,&
                                            &precipitation1,&
                                            &surface_temperature1,&
                                            &SNH, SMB, MELT, ACC, REFR, Albedo
  logical, dimension(:,:,:), allocatable :: mask1
  real, dimension(:,:), allocatable      :: snh_tmp, snh_lastyear
  real, dimension(:,:), allocatable      :: x, y
  real, dimension(:), allocatable        :: t
  integer :: xdim,ydim,tdim
  ! Define parameters
  real, parameter :: slim   = 7.       ! melt/refreeze threshold
  real, parameter :: Ans    = .85     ! Albedo for new snow
  real, parameter :: Ads    = .75      ! Albedo for dry snow
  real, parameter :: Aws    = .6     ! Albedo for wet snow
  real, parameter :: Ai     = .5     ! Albedo for ice
  real, parameter :: taucs  = .8
  real, parameter :: tauoc  = .3
  real, parameter :: epsacs = .7     ! emissivity of atm for sunny days
  real, parameter :: epsaoc = .95     ! emissivity of atm for overcasted days
  real, parameter :: epsi   = .95     ! emissivity of ice
  real, parameter :: beta   = 10.     ! turbulent heat transfer coeff
  real, parameter :: bolz   = 5.67e-8 ! Stefan-Boltzmann constant
  real, parameter :: T0     = 273.15  ! melt point in K
  real, parameter :: Tmin   = -6.5    ! background melt condition
  real, parameter :: pi     = 3.141592653
  real, dimension(12) :: mth = (/1,2,3,4,5,6,7,8,9,10,11,12/)
  real, dimension(12) :: mth_days = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  real, dimension(12) :: days = (/16,45,75,105,136,166,197,228,258,289,319,350/)
  integer, parameter :: mlen   = 12
  integer :: nlen
  real :: c1cs, c2cs, c1oc, c2oc
  real :: winkelns, winkelds, winkelws, winkeli

contains

  SUBROUTINE CALC_PARAMETER(c1cs, c2cs, c1oc, c2oc)
  ! ************************************************************************
  ! * This SUBTROUTINE calculates parameters c1 and c2                     *
  ! ************************************************************************
  real, intent(out) :: c1cs, c2cs, c1oc, c2oc
  ! c1 and c2 for clear sky conditions
  c1cs     = (epsi*epsacs*4*bolz*(T0**3)+beta)
  c2cs     = (-epsi+epsacs*epsi)*bolz*(T0**4)
  ! c1 and c2 for completely overcast conditions
  c1oc     = (epsi*epsaoc*4*bolz*(T0**3)+beta)
  c2oc     = (-epsi+epsaoc*epsi)*bolz*(T0**4)
  END SUBROUTINE CALC_PARAMETER

  SUBROUTINE VAR_ALLOCATE
  ! ************************************************************************
  ! * This SUBTROUTINE allocate memory for variables                       *
  ! ************************************************************************
    allocate (shortwave_radiation_downward(xlen,ylen,tlen))
    allocate (surface_temperature(xlen,ylen,tlen))
    allocate (precipitation(xlen,ylen,tlen))
    allocate (mask(xlen,ylen,tlen))
    allocate (mask1(xlen,ylen,tlen))
    allocate (snh_tmp(xlen,ylen))
    allocate (snh_lastyear(xlen,ylen))
    allocate (x(xlen,ylen))
    allocate (y(xlen,ylen))
    allocate (t(tlen))
    nlen = tlen/12
    allocate (lat(xlen,ylen,mlen,nlen))
    allocate (shortwave_radiation_downward1(xlen,ylen,mlen,nlen))
    allocate (surface_temperature1(xlen,ylen,mlen,nlen))
    allocate (precipitation1(xlen,ylen,mlen,nlen))
    allocate (SNH(xlen,ylen,mlen,nlen))
    allocate (SMB(xlen,ylen,mlen,nlen))
    allocate (MELT(xlen,ylen,mlen,nlen))
    allocate (ACC(xlen,ylen,mlen,nlen))
    allocate (REFR(xlen,ylen,mlen,nlen))
    allocate (Albedo(xlen,ylen,mlen,nlen))
  END SUBROUTINE VAR_ALLOCATE

  SUBROUTINE READ_NETCDF
    ! ************************************************************************
    ! * This SUBTROUTINE reads in atmosphere fields and ice-mask             *
    ! * note that this subroutine is orginally generated for coupling        *
    ! *     AWICM-PISM (developed by Paul)                                   *
    ! * for other input, units needs to be carefully checked                 *
    ! * atmosphere fields include :                                          *
    ! *   surface_temperature: units K                                       *
    ! *   shortwave_radiation_downward: units                                *
    ! *   precipitation: units kg m-2 second-1
    ! ************************************************************************
    IMPLICIT NONE

    integer :: status, ncid, varid
    integer :: i, j
    character (len = 20) :: precipitation_units,&
            &shortwave_radiation_downward_units,&
            &surface_temperature_units

    include 'netcdf.inc'

    status=nf_open(trim(filename_in), nf_nowrite, ncid)

    ! get shortwave_radiation_downward
    status=nf_inq_varid(ncid, shortwave_radiation_downward_varname, varid)
    status=nf_get_var_real(ncid, varid, shortwave_radiation_downward)
    shortwave_radiation_downward1 = reshape( shortwave_radiation_downward, (/xlen,ylen,mlen,nlen/) )
    ! write(*,*) "shortwave_radiation_downward1",shortwave_radiation_downward1(135,58,:,:)

    ! get surface_temperature
    status=nf_inq_varid(ncid, temperature_varname, varid)
    status=nf_get_var_real(ncid, varid, surface_temperature)
    status=nf_get_att_text(ncid, varid, "units", surface_temperature_units)
    !write(*,*) "surface_temperature_units is ",surface_temperature_units
    if (status.ne.nf_noerr) then
        write(*,*) 'WARNING: units of surface_temperature not found'
        write(*,*) 'Take K as default units'
        write(*,*) 'Possibly wrong result! Plz double check your result!'
        write(*,*) 'expected units for surface temperature: K or degC'
        surface_temperature = surface_temperature - 273.15
    else
        if (surface_temperature_units=='K') then
                !write(*,*) 'input surface_temperature units is K'
                surface_temperature = surface_temperature - 273.15
        elseif (surface_temperature_units=='degC') then
                !write(*,*) 'input surface_temperature units is degC'
        else
                write(*,*) 'WARNING: Invalid units of surface temperature'
                write(*,*) 'Invalid units',surface_temperature_units
                write(*,*) 'Expected units: K or degC'
                write(*,*) 'Here we will take K as default units'
                write(*,*) 'Possibly wrong result! Plz double check your result!'
                surface_temperature = surface_temperature - 273.15
        end if
    end if
    surface_temperature1 = reshape( surface_temperature, (/xlen,ylen,mlen,nlen/) )
    ! write(*,*) "surface_temperature1",surface_temperature1(135,58,:,:)

    ! get precipitation
    status=nf_inq_varid(ncid, precipitation_varname, varid)
    status=nf_get_var_real(ncid, varid, precipitation)
    status=nf_get_att_text(ncid, varid, "units", precipitation_units)
    !write(*,*) "precipitation_units is ",precipitation_units
    if (status.ne.nf_noerr) then
        write(*,*) 'WARNING: units of precipitation not found'
        write(*,*) 'Take kg m-2 second-1 as default units'
        write(*,*) 'Possibly wrong result! Plz double check your result!'
        write(*,*) 'expected units for precipitation: kg m-2 second-1,&
                &or mm day-1'
        precipitation = precipitation*24.*60.*60.
    else
        if (precipitation_units=='kg m-2 second-1') then
                !write(*,*) 'input precipitation units is kg m-2 second-1'
                precipitation = precipitation*24.*60.*60.
        elseif (precipitation_units=='mm day-1') then
                !write(*,*) 'input precipitation units is mm day-1'
        else
                write(*,*) 'WARNING: Invalid units of precipitation'
                write(*,*) 'Invalid units',precipitation_units
                write(*,*) 'Expected units: "kg m-2 second-1" or "mm day-1"'
                write(*,*) 'Here we will take "kg m-2 second-1" as default units'
                write(*,*) 'Possibly wrong result! Plz double check your result!'
                precipitation = precipitation*24.*60.*60.
        end if
    end if
    precipitation1 = reshape( precipitation, (/xlen,ylen,mlen,nlen/) )
    ! write(*,*) "precipitation",precipitation1(135,58,:,:)

    ! get ice mask
    !status=nf_inq_varid(ncid, mapping_varname, varid)
    !status=nf_get_var_real(ncid, varid, mask)
    !mask1 =  (mask(:,:,:) > 5)
    ! where (mask(:,:,:) > 100) mask1=1
    mask1=.True.

    ! get longitude
    status=nf_inq_varid(ncid, longitude_varname, varid)  !inquire var
    status=nf_get_var_real(ncid,varid,x)  ! to get x

    ! get latitude
    status=nf_inq_varid(ncid, latitude_varname, varid)  !inquire var
    status=nf_get_var_real(ncid,varid,y)
    ! spread latitude into four dimensions
    do i = 1, mlen
        do j = 1, nlen
           lat(:,:,i,j)=y
        end do
    end do

    !write(*,*) "latitude:",y

    ! get time
    status=nf_inq_varid(ncid, time_varname, varid)  !inquire var
    status=nf_get_var_real(ncid,varid,t)

    ! close
    status=nf_close(ncid)

    END SUBROUTINE


    SUBROUTINE READ_DATA

      ! write(*,*) "Reading from namelist..."
      CALL READ_NAMELIST

      ! write(*,*) "Allocating memories for variables..."
      CALL VAR_ALLOCATE

      ! write(*,*) "Reading in atmosphere fields...."
      CALL READ_NETCDF

    END SUBROUTINE READ_DATA

END MODULE MOD_DATA
