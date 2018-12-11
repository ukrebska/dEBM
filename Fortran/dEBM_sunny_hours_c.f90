subroutine dEBM_sunny_hours_c(mth,lat,elev,obl)
!******************************************************************************
! calculates the time, the sun is above a certain elevation angle
! elevation needs to be a scalar
! lat may be scalar, vector or matrix
! declination needs to be scalar, or same size as lat
! calls functions:
!                  dEBM_decl
!                  dEBM_sinphisind
!                  dEBM_cosphicosd(lat,decl)
!                  dEBM_sol_flux_fact
!                  dEBM_hangle2hours

USE DATA_MODULE
call calc_parameter

! for the middle of each month we calculate the number of days since
!  new year in a Julian calendar (very simple)
days=round(.5*(cumsum([31;28;31;30;31;30;31;31;30;31;30;31])+cumsum([0;31;28;31;30;31;30;31;31;30;31;30])))

k     = length(mth);
[m,n] = size(lat);
decl  = dEBM_decl(days(mth),obl);    ! mid-month decl of
decl  = reshape(decl,1,1,k);           ! just inflating declination
decl  = decl(ones(m,1),ones(n,1),:);   ! to mxnxk
lat   = lat(:,:,ones(k,1));            ! inflating lat to mxnxk
elev  = elev*pi/180;

!products of sin(lat)*sin(decl) and cos(lat)*cos(decl)  is needed more than once
sinphisind=sin(pi/180*lat).*sin(pi/180*decl)
cosphicosd=COS(pi/180*lat).*COS(pi/180*decl)

! hourangle of sun rise and factor which relates solar flux density and
! daily insolation SW
ha_0 = real(acos(-sinphisind/cosphicosd))!
sol_flux_fact_0=(sin_phi_sin_d*ha+cos_phi_cos_d*sin(ha_0))/pi!

! hourangle of sun rising above elev
ha_elev = real(acos((sin(elev)-sinphisind)/cosphicosd))!

! duration of melt period in hours
hours = ha_elev*24/pi;

! proportion of daily insolation which is effective during melt period
sol_flux_fact_e = (sinphisind*ha_elev+cosphicosd*sin(ha_elev))/pi!
q          = sol_flux_fact_e/sol_flux_fact_0!

end subroutine dEBM_sunny_hours_c
