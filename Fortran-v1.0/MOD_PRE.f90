MODULE MOD_PRE
! ************************************************************************
! * This MODULE reads info from namelist_pre_new                         *
! ************************************************************************

  integer       :: xlen, ylen, tlen
  real          :: stddev, obliquity
  character*300 :: filename_in
  character*20  :: precipitation_varname,&
                     &temperature_varname,&
                     &mapping_varname,&
                     &shortwave_radiation_downward_varname,&
                     &longitude_varname,&
                     &latitude_varname,&
                     &time_varname

contains

  SUBROUTINE READ_NAMELIST

    namelist /debm/ xlen, ylen,tlen,&
                    &filename_in,&
                    &precipitation_varname,&
                    &temperature_varname,&
                    &shortwave_radiation_downward_varname,&
                    &mapping_varname,&
                    &longitude_varname,latitude_varname,time_varname,&
                    &stddev, obliquity

    ! read information from namelist
    open(10, file="namelist.debm", status='old' )
    read(10, nml=debm)

   ! write(*,*) obliquity

    close(10)

 END SUBROUTINE READ_NAMELIST

END MODULE MOD_PRE
