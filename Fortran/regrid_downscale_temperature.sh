#!/bin/sh

INPUT_FILE_ATMOSPHERE="/Users/xushan/Downloads/dEBM-master/pi477_u_MM_2601_270001.01_echam.nc"
ATMOSPHERE_NAME_TEMPERATURE="temp2"
ATMOSPHERE_NAME_ELEVATION="gld"
ATMOSPHERE_NAME_PRECIPITATION="aprl"
INPUT_FILE_PISM="/Users/xushan/Code/dEBM/Fortran/data/pism.nc"
DOWNSCALING_LAPSE_RATE=-0.007        # FIXME: should this be 7, or 7/1000?
GAMMA=0.07

# ******************************************************************************
# getting high resolution grid information
cdo griddes  \
        ${INPUT_FILE_PISM} > \
        grid.txt

# regridding temperature &  bilinearly from ${INPUT_FILE_ATMOSPHERE}"
cdo -s -f nc  \
        -remapbil,grid.txt   \
        -selvar,${ATMOSPHERE_NAME_TEMPERATURE} \
        ${INPUT_FILE_ATMOSPHERE} \
        temperature_lo.nc

# regridding low resolution topo from ${INPUT_FILE_ATMOSPHERE} into high resolution"
cdo -s -f nc  \
        -remapbil,grid.txt   \
        -selvar,${ATMOSPHERE_NAME_ELEVATION} \
        ${INPUT_FILE_ATMOSPHERE} \
        elev_lo.nc

# generating elevation difference
INPUT_FILE_NAME_BASE_pism=$(basename $INPUT_FILE_PISM)
if cdo -s pardes ${INPUT_FILE_PISM} | grep -q usurf; then
    cdo -s \
        -seltimestep,-1 \
        -selvar,usurf ${INPUT_FILE_PISM} \
        "${INPUT_FILE_NAME_BASE_pism%.*}"_usurf.nc
    elev_hi_file="${INPUT_FILE_NAME_BASE_pism%.*}"_usurf.nc
else
    cdo -s -add \
        -selvar,thk ${INPUT_FILE_PISM} \
                    -setrtoc,-10000,0,0 \
        -selvar,topg ${INPUT_FILE_PISM} \
        "${INPUT_FILE_NAME_BASE_pism%.*}"_thk_plus_topg.nc
    elev_hi_file="${INPUT_FILE_NAME_BASE_pism%.*}"_thk_plus_topg.nc
fi
mv $elev_hi_file elev_hi.nc
cdo -s sub elev_hi.nc elev_lo.nc elev_hi_minus_lo.nc

# downscaling temperature
cdo -s -f nc \
        -add temperature_lo.nc \
        -mulc,${DOWNSCALING_LAPSE_RATE} \
        elev_hi_minus_lo.nc \
        temperature_hi.nc
mv temperature_hi.nc temperature.nc

# ******************************************************************************
# Now we do precipitation
# The following is after Quiquet et al, The Cryosphere, 2012
# https://www.the-cryosphere.net/6/999/2012/
# The factor gamma is chosen after Huybrechts 2002, but is poorly constrained (see Charbit et al. 2002)

# regridding precipitation &  bilinearly from ${INPUT_FILE_ATMOSPHERE}"
cdo -s -f nc  \
        -remapbil,grid.txt   \
        -selvar,${ATMOSPHERE_NAME_PRECIPITATION} \
        ${INPUT_FILE_ATMOSPHERE} \
        precipitation_lo.nc

# downscaling precipitation
cdo -s \
        -mul precipitation_lo.nc \
        -exp \
          -mulc,${GAMMA} \
          -sub temperature.nc temperature_lo.nc \
        precip_hi.nc
mv precip_hi.nc precipitation.nc

# make clean
rm grid.txt
rm temperature_lo.nc elev_lo.nc elev_hi.nc elev_hi_minus_lo.nc
rm precipitation_lo.nc
