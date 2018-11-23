# dEBM
dEBM is a surface melt scheme to couple ice and climate models in paleo applications. 


`dEBMmain.m` takes near surface air temperature (˚C) and downward short wave radiation (monthly mean of daily means flux in W/m2), albedo and latitude as input, and calls all the necessary routines to calculate melt rates in mm/day.


Note that (for simplicity) the calendar is modern and declination for a given date is only approximated in a very simple way 
(it does not matter too much, but surely this can be improved).

A full surface mass balance model is on its way and may be published soon.

Citation: 
> Krebs-Kanzow, U., Gierz, P., and Lohmann, G., Brief communication: An Ice surface melt scheme including the diurnal cycle of solar radiation, The Cryosphere Discuss., https://doi.org/10.5194/tc-2018-130, accepted for publication

See also `CITATION.md` for a BibTeX entry

(C) Uta Krebs-Kanzow, Alfred Wegener Institute, Bremerhaven, Germany, 2018
