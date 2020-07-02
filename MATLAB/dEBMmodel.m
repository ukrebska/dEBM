function melt=dEBMmodel(A, swd, PDD,hours,q,c1,c2,MC)
% please cite:
% Krebs-Kanzow, U., Gierz, P., and Lohmann, G.:
% Brief communication: An Ice surface melt scheme including the diurnal cycle
% of solar radiation,
% The Cryosphere Discuss.,
% https://doi.org/10.5194/tc-2018-130, in press
% (C) Uta Krebs-Kanzow, Alfred Wegener Institute, Bremerhaven, Germany, 2018
%*****************************************************************************
% A is albedo
% swd is daily short wave insolation at the surface
% PDD is an approx. of melt-period temperature (e.g. PDD(T,3.5))
% hours is the duration of the daily melt period in hours
% q is the proportion of sw which is effective during the melt period
% MC is (logic) the background melt condition (e.g. T>-6.5)  
% melt is the melt rate (mm/day) 
% Lf is latent heat of fusion
Lf=3.34e5;
melt=MC.*max(0,((1-A).*swd.*24.*q)+(c2+c1*PDD).*hours)/Lf*60*60;
