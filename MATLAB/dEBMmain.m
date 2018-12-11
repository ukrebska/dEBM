function melt=dEBMmain(lat,swd,T,A)
% please cite:
% Krebs-Kanzow, U., Gierz, P., and Lohmann, G.:
% Brief communication: An Ice surface melt scheme including the diurnal cycle
% of solar radiation,
% The Cryosphere Discuss.,
% https://doi.org/10.5194/tc-2018-130, in press
% (C) Uta Krebs-Kanzow, Alfred Wegener Institute, Bremerhaven, Germany, 2018
%******************************************************************************
% A,lat,swd,PDD are here expexted to be 
% fields of size (xdim x ydim x 12 x number of years) 
% 

% parameters defined for dEBM 
[xdim,ydim,bla,NY] = size(swd);
dEBMmelt           = zeros(xdim,ydim,bla,NY); (initialize)
epsa   = .76;     % atm. emissivity
epsi   = .95;     % emissivity of ice
beta   = 10;      % turbulent heat transfer coeff
bolz   = 5.67e-8; % Stefan-Boltzmann constant
T0     = 273.15;  % melt point in K
Tmin   = -6.5;    % background melt condition
c1     = (epsa*4*bolz*T0^3+beta);   
c2     = (-epsi+epsa*epsi)*bolz*T0^4;
elev   = 17.5
obl    =23.44;
mth    = 1:12     

% generate parameters that depend on latitude and month 
[hours,q,bla] = dEBM_sunny_hours_c(mth,lat,elev,obl);

% inflate parameters to create (xdim x ydim x 12 x number of years) 
hours = hours(:,:,:,ones(NY,1));       
q     = q(:,:,:,ones(NY,1));

% calculte PDD 
PDD   = PDD4(T,3.5);

% melt condition: ice is warm enough to warm to melting point during daytime	
% MC    = (T(:,:,mth,i)>Tmin);    %original
MC    = (T(:,:,mth,:)>Tmin);      % modified by Shan  

% calculate monthly melt rates		
for i=1:NY
   for mth=1:12	
       dEBMmelt(:,:,mth,i) = dEBMmodel(A(:,:,mth,i),...
				       swd(:,:,mth,i),...
				       PDD(:,:,mth,i),...
				       hours(:,:,mth,i),...
                                       q(:,:,mth,i),c1,c2,MC(:,:,mth,i));     % modified by Shan
% 				       q(:,:,mth,i),c1,c2,MC);                %original				       
   end
end
