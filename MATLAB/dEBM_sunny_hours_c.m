function [mp_hours,q_elev,sol_flux_fact_0]=dEBM_sunny_hours_c(mth,lat,elev,obl)
% please cite:
% Krebs-Kanzow, U., Gierz, P., and Lohmann, G.:
% Brief communication: An Ice surface melt scheme including the diurnal cycle
% of solar radiation,
% The Cryosphere Discuss.,
% https://doi.org/10.5194/tc-2018-130, in press
% (C) Uta Krebs-Kanzow, Alfred Wegener Institute, Bremerhaven, Germany, 2018
%******************************************************************************
% calculates the time, the sun is above a certain elevation angle
% elevation needs to be a scalar
% lat may be scalar, vector or matrix
% declination needs to be scalar, or same size as lat
% calls functions: dEBM_day
%                  dEBM_decl
%                  dEBM_sinphisind
%                  dEBM_cosphicosd(lat,decl)
%                  dEBM_sol_flux_fact
%                  dEBM_hangle2hours   
days  = dEBM_day();   % day of the year at mid-month for each month of the year
k     = length(mth);
[m,n] = size(lat); 
decl  = dEBM_decl(days(mth),obl);    % mid-month decl of 
decl  = reshape(decl,1,1,k);           % just inflating declination       
decl  = decl(ones(m,1),ones(n,1),:);   % to mxnxk
lat   = lat(:,:,ones(k,1));            % inflating lat to mxnxk
elev  = elev*pi/180;

%products of sin(lat)*sin(decl) and cos(lat)*cos(decl)  is needed more than once
sinphisind=dEBM_sinphisind(lat,decl);
cosphicosd=dEBM_cosphicosd(lat,decl);

% hourangle of sun rise and factor which relates solar flux density and
% daily insolation SW
ha_0 = real(acos(-sinphisind./cosphicosd));
sol_flux_fact_0 = dEBM_sol_flux_fact(sinphisind,cosphicosd,ha_0);

% hourangle of sun rising above elev  
ha_elev = real(acos((sin(elev)-sinphisind)...
		            ./cosphicosd));
% duration of melt period in hours	    

mp_hours = dEBM_hangle2hours(ha_elev);       

% proportion of daily insolation which is effective during melt period
sol_flux_fact_e = dEBM_sol_flux_fact(sinphisind,cosphicosd,ha_elev);
q_elev          = sol_flux_fact_e./sol_flux_fact_0;

