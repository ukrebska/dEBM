function sinphisind=dEBM_sin_phi_sin_d(lat,decl)
% just here to avoid unnecessary re-computation   
sinphisind=sin(pi/180*lat).*sin(pi/180*decl);
