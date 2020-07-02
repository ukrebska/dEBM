########################################
# Set up how to compile the source files
########################################


# Add the source files
SET(DEBM_src ${SRCDEBM}/dEBMmain.f90
             ${SRCDEBM}/MOD_DATA.f90
             ${SRCDEBM}/MOD_MAIN.f90
             ${SRCDEBM}/MOD_OUTPUT.f90
             ${SRCDEBM}/MOD_PRE.f90
             ${SRCDEBM}/MOD_RESTART.f90
)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${DEBMEXE} ${DEBM_src})

# Allow NetCDF:
TARGET_INCLUDE_DIRECTORIES(${DEBMEXE} PRIVATE ${NETCDF_Fortran_INCLUDE_DIRECTORIES})
TARGET_LINK_LIBRARIES(${DEBMEXE} PRIVATE ${NETCDF_Fortran_LIBRARIES})

# PG: Probably wrong:
SET(CMAKE_INSTALL_PREFIX /usr/local)
INSTALL(TARGETS ${DEBMEXE} RUNTIME DESTINATION bin)
