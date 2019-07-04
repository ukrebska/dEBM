PROGRAM MAIN
  USE MOD_PRE
  USE MOD_DATA
  USE MOD_MAIN
  USE MOD_OUTPUT

  IMPLICIT NONE

  integer :: n

  write(*,*) "********* Begin *********"

  write(*,*) "Starting calculating dEBM..."
  CALL READ_DATA

  snh_tmp   = 0.
  snh_lastyear = 0.
  SNH      = 0.
  SMB      = 0.
  MELT     = 0.
  ACC      = 0.
  REFR     = 0.
  Albedo   = 0.

  ! write(*,*) "calculates the first year...."
  CALL dEBM_core(surface_temperature1(:,:,:,1), shortwave_radiation_downward1(:,:,:,1), precipitation1(:,:,:,1), snh_tmp, snh_lastyear, lat(:,:,:,1), mask1, obliquity, SNH(:,:,:,1), SMB(:,:,:,1), MELT(:,:,:,1), ACC(:,:,:,1), REFR(:,:,:,1), Albedo(:,:,:,1))
  snh_tmp      = SNH(:,:,12,1)
  snh_lastyear = SNH(:,:,9,1)
!
  ! write(*,*) "SNH",SNH(376,240,:,1)
  ! write(*,*) "SMB",SMB(376,240,:,1)
  ! write(*,*) "MELT",MELT(376,240,:,1)
  ! write(*,*) "ACC",ACC(376,240,:,1)
  ! write(*,*) "REFR",REFR(376,240,:,1)
  ! write(*,*) "Albedo",Albedo(376,240,:,1)

  do n = 1, nlen
      write(*,*) "calculates year",n
      CALL dEBM_core(surface_temperature1(:,:,:,n), shortwave_radiation_downward1(:,:,:,n), precipitation1(:,:,:,n), snh_tmp, snh_lastyear, lat(:,:,:,n), mask1, obliquity, SNH(:,:,:,n), SMB(:,:,:,n), MELT(:,:,:,n), ACC(:,:,:,n), REFR(:,:,:,n), Albedo(:,:,:,n))
      snh_tmp      = SNH(:,:,12,n)
      snh_lastyear = SNH(:,:,9,n)
  end do

  ! write(*,*) "SNH",SNH(376,240,:,1)
  ! write(*,*) "SMB",SMB(376,240,:,1)
  ! write(*,*) "MELT",MELT(376,240,:,1)
  ! write(*,*) "ACC",ACC(376,240,:,1)
  ! write(*,*) "REFR",REFR(376,240,:,1)
  ! write(*,*) "Albedo",Albedo(376,240,:,1)

  !write(*,*) "latitude output is ",y
  CALL OUTPUT_NETCDF(x, y, SNH, SMB, MELT, ACC, REFR, Albedo)

  write(*,*) "********* Done. *********"

END PROGRAM MAIN
