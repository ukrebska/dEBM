MODULE MOD_PRE
  integer       :: xlen, ylen,tlen
  character*300 :: filename_in
  character*20  :: albedo_varname,temperature_varname,&
                     &shortwave_radiation_downward_varname,longitude_varname,latitude_varname,time_varname
  real          :: stddev
contains
  SUBROUTINE READ_NAMELIST
    namelist /na/ xlen, ylen,tlen,&
                    &filename_in,albedo_varname,temperature_varname,shortwave_radiation_downward_varname,&
                    &longitude_varname,latitude_varname,time_varname,stddev
    open(10, file="./namelist_pre_new.txt", status='old' )
    read (10, nml=na)
    close (10)
 END SUBROUTINE READ_NAMELIST
END MODULE MOD_PRE
