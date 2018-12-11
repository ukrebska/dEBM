MODULE DATA_MODULE
  IMPLICIT NONE
! input variables
real, dimension(:,:,:,:), allocatable :: swd
integer :: xdim, ydim, tdim, NY
! parameters defined for dEBM
real, parameter :: epsa   = .76     ! atm. emissivity
real, parameter :: epsi   = .95     ! emissivity of ice
real, parameter :: beta   = 10      ! turbulent heat transfer coeff
real, parameter :: bolz   = 5.67e-8 ! Stefan-Boltzmann constant
real, parameter :: T0     = 273.15  ! melt point in K
real, parameter :: Tmin   = -6.5    ! background melt condition
real :: c1, c2
real, parameter :: elev   = 17.5
real, parameter :: obl    = 23.44
integer, dimension(12) :: mth = (/1,2,3,4,5,6,7,8,9,10,11,12/)
integer, dimension(12) :: days = (/16,45,75,105,136,166,197,228,258,289,319,350/)
real, paramter :: Lf=3.34e5         ! Lf is latent heat of fusion
END MODULE DATA_MODULE



SUBROUTINE CALC_PARAMETER
  USE DATA_MODULE
  !
  c1     = (epsa*4*bolz*(T0**3)+beta)
  c2     = (-epsi+epsa*epsi)*bolz*(T0**4)
  ! allocate memory
  xdim   = size(swd,dim=1)
  ydim   = size(swd,dim=2)
  tdim   = size(swd,dim=3)
  NY     = size(swd,dim=4)
  allocate (swd(xdim,ydim,tdim,NY))
  !
  real, dimension(xdim,ydim) :: lat
  real, dimension(xdim,ydim,tdim,NY) :: A
  real, dimension(xdim,ydim,tdim,NY) :: PDD
  real, dimension(xdim,ydim,tdim,NY) :: dEBMmelt
  real, dimension(xdim,ydim,tdim,NY) :: hours
  real, dimension(xdim,ydim,tdim,NY) :: q
  real, dimension(xdim,ydim,tdim,NY) :: MC=0
  RETURN
END SUBROUTINE CALC_PARAMETER
