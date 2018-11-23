function decl = dEBM_decl(days,obliqu)
% so far this is a very simple calculation, which assumes eccentricity
%  is 1... and spring equinox 
% is on March 20th.  
% Values get more inaccurate towards autumn, lap years are neglected, too
% In the present day orbital configuration, the Earth orbits slower in boreal summer and decl=23.44*sin(pi/180*186/180*(days(:)-79)) is 
% more accurate for the northern hemisphere but not a good choice for the southern hemisphere...
% uncomment the following to see the difference
% figure
% hold
% plot(1:365,-23.44*sin(pi/180*180/186*((1:365)-79)),'k') % correct number of days between equinoxes
% plot(1:365,-23.44*sin(pi/180*360/365*((1:365)-79)),'r') % eccentricity = 1
% plot([79,79],[-24,24]) % spring equinox
% plot([266,266],[-24,24]) % autumn equinox
% plot([172,172],[-24,24]) % autumn equinox
decl=obliqu*sin(pi/180*360/365*(days(:)-79));   % 79 is the number of days since New Year for the spring equinox

 
