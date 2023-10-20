#!/bin/bash
#                
##module load gcc/11.2.0-gcc-11.2.0 #gcc/6.4.0
#module load intel-oneapi-compilers
module load nco/5.0.6-gcc-11.2.0 intel-oneapi-compilers/2022.0.1-gcc-11.2.0 intel-oneapi-mkl/2022.0.1-gcc-11.2.0 openmpi/4.1.2-intel-2021.5.0 netcdf-fortran/4.5.3-openmpi-4.1.2-intel-2021.5.0
##module load netcdf-fortran/4.5.3-gcc-11.2.0
#fortran_compiler=ifort #gfortran
##nc_config=/sw/spack-levante.2022-02-17/netcdf-fortran-4.5.3-jlxcfz/bin/nf-config #/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.4-gcc64/bin/nf-config
#nc_config=nc-config #/sw/spack-levante/netcdf-fortran-4.5.3-jlxcfz//bin/nf-config
nc_config=/sw/spack-levante/netcdf-fortran-4.5.3-pywf2l/bin/nf-config
export NETCDF_LIB=$($nc_config --flibs)
#export NETCDF_LIB=$($nc_config --libdir)
export NETCDF_INC=$($nc_config --includedir)

NETCDFFROOT=/sw/spack-levante/netcdf-fortran-4.5.3-pywf2l/
export PATH=${NETCDFFROOT}:${NETCDFFROOT}/include:${NETCDFFROOT}/lib:$PATH
export LD_LIBRARY_PATH=${NETCDFFROOT}/lib:${NETCDFROOT}/lib:$LD_LIBRARY_PATH

thisdir=$(dirname $0)
thisdir=$(readlink -f $thisdir)
old_dir=$(pwd)

cd ${thisdir}/Fortran/

make

rm -rf ${thisdir}/_build/
mkdir -p ${thisdir}/_build/



#echo $fortran_compiler \
#        ${thisdir}/calnoro.f90 \
#        ${thisdir}/grid_noro.f90 \
#        -o ${thisdir}/_build/calnoro \
#        $NETCDF_LIB $NETCDF_INCLUDE 
#$fortran_compiler \
#        MOD_PRE.f90 MOD_DATA.f90 MOD_MAIN.f90 MOD_OUTPUT.f90 dEBMmain.f90 \
#        -o ${thisdir}/_build/dEBMmain \
#        $NETCDF_LIB $NETCDF_INCLUDE 
