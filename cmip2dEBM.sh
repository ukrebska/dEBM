#!/bin/bash
weightfile=./Weights.AWICM_TO_5km.bil.nc
griddes=5km_grid.txt
idir=.
ifile=Climatology_AWIESM_1959_2000.nc

cdo -f nc -t echam6 -genbil,$griddes ${ifile} $weightfile
echo "${year}!"
ofile=test_inpute_dEBM_5km.nc
#cdo -merge ${idir}/tas_${fbase}_${year}01-${year}12.nc ${idir}/rsds_${fbase}_${year}01-${year}12.nc ${idir}/rlds_${fbase}_${year}01-${year}12.nc clt/clt_${fbase}_${year}01-${year}12.nc ${idir}/pr_${fbase}_${year}01-${year}12.nc ${idir}/rsdt_${fbase}_${year}01-${year}12.nc AWI_ESM_geosp.nc merge_input.nc 
cdo -s -O -f nc -t echam6 -copy -expr,'srad0d=rsdt; emiss=rlds/(5.6703744e-8*tas^4); swd=rsds; lwd=rlds; temp2=tas; precip=pr*86400; aclcov=clt; orog=geosp/9.81;' ${ifile} tmp.nc
cdo -s -remap,$griddes,$weightfile tmp.nc ${ofile}
rm tmp.nc
