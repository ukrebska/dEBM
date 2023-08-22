#!/bin/bash

#SBATCH --job-name=smb   # Specify job name
#SBATCH --output=smb.o%j    # name for standard output log file
#SBATCH --error=smb.e%j     # name for standard error output log file
#SBATCH --partition=shared      # Specify partition name
#SBATCH --ntasks=36             # Specify max. number of tasks to be invoked
#SBATCH --time=08:00:00        # Set a limit on the total run time
#SBATCH --account=ba0989      # Charge resources on this project account
#SBATCH --mem=6GB



set -e


module load nco
# default
cdo=$(which cdo)
python3=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin/python3

#### for dEBM, calnoro jsbach init binary
#NETCDFFROOT=/sw/spack-levante/netcdf-fortran-4.5.3-k6xq5g/
NETCDFFROOT=/sw/spack-levante/netcdf-fortran-4.5.3-pywf2l/
export PATH=${NETCDFFROOT}:${NETCDFFROOT}/include:${NETCDFFROOT}/lib:$PATH
export LD_LIBRARY_PATH=${NETCDFFROOT}/lib:${NETCDFROOT}/lib:$LD_LIBRARY_PATH




ECHAM_TO_ISM_multiyear_mean=1  #True
FESOM_TO_ISM_multiyear_mean=1  # True
number_of_years_for_forcing=3  # usually consistent with chunk size,


#dEBM_binary=/home/a/a270075/model_codes/pism/dEBM/Fortran-v1.0/dEBMmain
dEBM_binary=/home/a/a270124/model_codes/dEBM/Fortran/dEBMmain

PISM_griddes=/work/ba1066/a270124/pool_pism/grids/greenland/pismr_greenland_5km.griddes
RES_pism=5km

SMB_scheme='dEBM'
DOWNSCALE_scheme='bil'
LAPSERATE_T='linear' # 'None'
LAPSERATE_T_value=0  #-5
LAPSERATE_P='None'  # 'desertification'


# OCEAN_scheme='pico'
# PISM_ocean_basin=/home/a/a270075/ba0989/pool/pism/basins/Basins_NH_20km_for_PICO.nc

ATMOS_startJuly=False


# METHOD_total_ice_mass_loss_flux=2  #3
# HOSING_CORRECTION=0
##############33
FUNCTION_PATH=/home/a/a270124/esm_tools/runscripts/awicm-pism/coupling
COUPLE_DIR=$(pwd)
WORK_DIR=$(pwd)
source ${FUNCTION_PATH}/coupling_echam2ice.functions
source ${FUNCTION_PATH}/coupling_atmosphere2pism.functions


EXP_ID=tran20-38PISM   # define
file_start_topg=/home/a/a270075/ba0989/awiesm_pism/pool_pism/pism_topg_thk_heatflx_base21ice38pism.nc  #/home/a/a270075/ba0989/awiesm_pism/pool_pism/pism_topg_thk_bheatflx_glac1d_38ka.nc  #/home/a/a270075/ba0989/awiesm_pism/pool_pism/pism_topg_thk_bheatflx_37_5ka.nc  #define


expdir=/home/a/a270075/ba0989/awiesm_pism/experiments/tranx20/${EXP_ID}
DATA_DIR_echam=${expdir}/outdata/echam/



job_start=2000
job_end=2900


run_dy=5
run_end=$job_start
pism_end=99999

while [  $run_end -lt $job_end ]
do
    echo
  run_start=$run_end
  run_end=$[job_end<run_start+run_dy?job_end:run_start+run_dy]


  # define the pism restart file:
  if [ $run_start -eq $job_start ]; then
   INPUT_FILE_pism=$file_start_topg
  else
    INPUT_FILE_pism=${expdir}/restart/pism/${EXP_ID}_pismr_restart_${pism_start}0101-${pism_end}1231.nc
  fi
  ls $INPUT_FILE_pism

  pism_start=$(( pism_end + 1  ))
  pism_end=$(( pism_start + 99  ))

  # for the climate model year and the next pism model year
  # echo ${run_start}  ${run_end}
  # echo ${pism_start}  ${pism_end}

  CHUNK_START_DATE_echam=${run_start}-01-01T00:00:00
  CHUNK_END_DATE_echam=$(( run_end-1 ))-12-31T00:00:00
  CHUNK_START_DATE_pism=${pism_start}0101
  CHUNK_END_DATE_pism=${pism_end}1231

  echo  ${CHUNK_START_DATE_echam}  ${CHUNK_END_DATE_echam}
  echo ${CHUNK_START_DATE_pism}  ${CHUNK_END_DATE_pism}

  ## main calculation:
  namelist_echam=${expdir}/config/echam/namelist.echam_${CHUNK_END_DATE_echam:0:4}0101-${CHUNK_END_DATE_echam:0:4}1231
  echam2ice
  atmosphere2pism

done
