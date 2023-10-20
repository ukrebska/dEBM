dEBMmain:
	ifort -I$(NETCDF_INC) $(NETCDF_LIB) -lnetcdf -lnetcdff MOD_PRE.f90 MOD_DATA.f90 MOD_MAIN.f90 MOD_OUTPUT.f90 dEBMmain.f90 -o dEBMmain

