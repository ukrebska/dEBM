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
  real, parameter :: Lf     = 3.34e5         ! Lf is latent heat of fusion
  real, parameter :: pi     = 3.1415
END MODULE DATA_MODULE


FUNCTION dEBMmain(lat,swd,T,A)
  IMPLICIT NONE
  USE DATA_MODULE
  CALL CALC_PARAMETER
  ! generate parameters that depEND on latitude and month
  ![hours,q,bla] = dEBM_sunny_hours_c(mth,lat,elev,obl);
  CALL dEBM_sunny_hours_c(mth,lat,elev,obl)
  ! calculte PDD
  CALL PDD(T,3.5)
  ! melt condition: ice is warm enough to warm to melting point during daytime
  ! MC    = (T(:,:,mth,i)>Tmin);
  where (T>Tmin) MC = 1.
  ! calculate monthly melt rates
  ! dEBMmodel
  swd0 = (1-A)*swd*24*q!
  where (swd0 < 0) swd0 = 0.
  dEBMmelt=MC*swd0+(c2+c1*PDD)*hours)/Lf*60*60!
END FUNCTION dEBMmain

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


SUBROUTINE dEBM_sunny_hours_c(mth,lat,elev,obl)
  !******************************************************************************
  ! calculates the time, the sun is above a certain elevation angle
  ! elevation needs to be a scalar
  ! lat may be scalar, vector or matrix
  ! declination needs to be scalar, or same size as lat
  USE DATA_MODULE
  call calc_parameter
  real, dimension(xdim,ydim,tdim) :: decl
  real, dimension(xdim,ydim,tdim) :: lat_new
  real, dimension(xdim,ydim,tdim) :: sinphisind,cosphicosd
  real, dimension(xdim,ydim,tdim) :: ha_0,ha_elev,mp_hours
  real, dimension(xdim,ydim,tdim) :: sol_flux_fact_0,sol_flux_fact_e,q_elev
  ! for the middle of each month we calculate the number of days since
  !  new year in a Julian calENDar (very simple)
  !days=round(.5*(cumsum([31;28;31;30;31;30;31;31;30;31;30;31])+cumsum([0;31;28;31;30;31;30;31;31;30;31;30])))
  !integer, dimension(12) :: days = (/16,45,75,105,136,166,197,228,258,289,319,350/)
  !decl  = dEBM_decl(days(mth),obl);    ! mid-month decl of
  do i=1, xdim
    do j=1, ydim
     decl(i,j,:)   = obl*sin(pi/180*360/365*(days-79))!  ! 79 is the number of days since New Year for the spring equinox
    END do
  END do
  ! inflating lat to mxnxk
  do i=1,mth
    lat_new(:,:,i) = lat
  END do
  ! elev
  elev             = elev*pi/180!
  !products of sin(lat)*sin(decl) and cos(lat)*cos(decl)  is needed more than once
  sinphisind       = sin(pi/180*lat_new)*sin(pi/180*decl)!
  cosphicosd       = cos(pi/180*lat_new)*cos(pi/180*decl)!
  ! hourangle of sun rise and factor which relates solar flux density and
  ! daily insolation SW
  ha_0             = real(acos(-sinphisind/cosphicosd))!
  sol_flux_fact_0  = (sin_phi_sin_d*ha+cos_phi_cos_d*sin(ha_0))/pi!
  ! hourangle of sun rising above elev
  ha_elev          = real(acos((sin(elev)-sinphisind)/cosphicosd))!
  ! duration of melt period in hours
  mp_hours         = ha_elev*24/pi!
  do t=1,tdim
    hours(:,:,:,t) = mp_hours
  end do
  ! proportion of daily insolation which is effective during melt period
  sol_flux_fact_e  = (sinphisind*ha_elev+cosphicosd*sin(ha_elev))/pi!
  q_elev           = sol_flux_fact_e/sol_flux_fact_0!
  do t=1,tdim
    q(:,:,:,t) = q_elev
  end do
  RETURN
END SUBROUTINE dEBM_sunny_hours_c

SUBROUTINE PDD(T,stddev)
  ! PDD is approximated as described in
  ! Calov R and Greve R (2005)
  ! Correspondence. A semi-analytical solution for the positive degree-day model
  ! with stochastic temperature variations.
  ! J. Glaciol., 51 (172), 173–175
  ! (doi:10.3189/172756505781829601)
  !****************************************************************************
  ! T in ^oC!!!
  ! stddev can be spatially variable (same grid as T) or a constant.
  !T=T-273.15;
  integer :: stddev
  PDD=stddev/sqrt(2*pi)*exp(-(T**2)/(2*(stddev**2)))+T/2.*erfc(-T/(sqrt(2)*stddev))!
END SUBROUTINE PDD
