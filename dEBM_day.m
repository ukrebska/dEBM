function days = dEBM_day()
% for the middle of each month we calculate the number of days since
%  new year in a Julian calendar (very simple)
days=round(.5*(cumsum([31;28;31;30;31;30;31;31;30;31;30;31])+cumsum([0;31;28;31;30;31;30;31;31;30;31;30])));
