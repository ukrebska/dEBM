MODULE DATA_MODULE
  IMPLICIT NONE
  ! input variables
  real, dimension(:,:,:,:), allocatable :: swd,A,PDD,dEBMmelt,hours,q,MC
  real, dimension(:,:),     allocatable :: lat
  integer :: xdim, ydim, tdim, nyear
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
contains
    SUBROUTINE CALC_PARAMETER
      real :: c1, c2
      c1     = (epsa*4*bolz*(T0**3)+beta)
      c2     = (-epsi+epsa*epsi)*bolz*(T0**4)
      RETURN
    END SUBROUTINE CALC_PARAMETER
END MODULE DATA_MODULE

PROGRAM TEST
  USE DATA_MODULE
  !implicit none
  integer :: ierr,ncid,varid,varidx,varidy,len_file,err
  ! real  :: x(LLm,1)
  ! real  :: y(MMm,1)
  ! real  :: lon(LLm,MMm),lat(LLm,MMm)
  ! real  :: swd(LLm,MMm,TTm),A(LLm,MMm,TTm),T(LLm,MMm,TTmi)
  ! real  :: dEBMmelt(LLm,MMm,TTm)
  !
  character*300 :: file_in
  include 'netcdf.inc'
  file_in='/Users/xushan/Downloads/dEBM-master/pi477_u_MM_2601_270001.01_echam.nc'
  len_file=len_trim(file_in)
  ierr=nf_open(trim(file_in),nf_nowrite,ncid)
  ierr=nf_inq_varid(ncid,'rsuscs_na',varid)
  ierr=nf_get_var_real(ncid,varid,swd)

  xdim   = size(swd,dim=1)
  ydim   = size(swd,dim=2)
  tdim   = size(swd,dim=3)
  NY     = size(swd,dim=4)
  ! allocate (swd(xdim,ydim,tdim,NY))
  allocate (lat(xdim,ydim))
  ! allocate (PDD(xdim,ydim,tdim,NY))
  ! allocate (dEBMmelt(xdim,ydim,tdim,NY))
  ! allocate (hours(xdim,ydim,tdim,NY))
  ! allocate (q(xdim,ydim,tdim,NY))
  ! allocate (MC(xdim,ydim,tdim,NY))

  ! ierr=nf_inq_varid(ncid,'temp2',varid)
  ! ierr=nf_get_var_real(ncid,varid,T)
  ! ierr=nf_inq_varid(ncid,'alsoi',varid)
  ! ierr=nf_get_var_real(ncid,varid,A)
  ! ierr=nf_inq_varid(ncid,'lon',varidx)
  ! ierr=nf_inq_varid(ncid,'lat',varidy)
  ! ierr=nf_get_var_real(ncid,varidx,x)
  ! ierr=nf_get_var_real(ncid,varidy,y)
  ! CALL meshgrid(x,y,lon,lat)
!dEBMmelt=dEBMmain(lat,swd,T,A)
  !print *,swd(1,1,1),swd(154,4,1)
  !print *,lon
  !print *,lat
  ! CALL ALLOCATE_MEMORY()
  err=NF_CLOSE(ncid)
  write(*,*) xdim,ydim,tdim,NY
    !implicit none
    !
    !
END PROGRAM TEST

! SUBROUTINE ALLOCATE_MEMORY
!   USE DATA_MODULE
!   ! allocate memory
!   xdim   = size(swd,dim=1)
!   ydim   = size(swd,dim=2)
!   tdim   = size(swd,dim=3)
!   NY     = size(swd,dim=4)
!   allocate (swd(xdim,ydim,tdim,NY))
!   allocate (lat(xdim,ydim))
!   allocate (PDD(xdim,ydim,tdim,NY))
!   allocate (dEBMmelt(xdim,ydim,tdim,NY))
!   allocate (hours(xdim,ydim,tdim,NY))
!   allocate (q(xdim,ydim,tdim,NY))
!   allocate (MC(xdim,ydim,tdim,NY))
!   RETURN
! END SUBROUTINE ALLOCATE_MEMORY

! SUBROUTINE meshgrid(x, y, x2, y2)
! real(dp), intent(in) :: x(:), y(:)
! real(dp), intent(out) :: x2(:, :), y2(:, :)
! x2 = spread(x, 1, size(y))
! y2 = spread(y, 2, size(x))
! END SUBROUTINE
