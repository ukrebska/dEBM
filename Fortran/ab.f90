PROGRAM TEST
  implicit none
  include 'netcdf.inc'
  integer :: ierr,ncid,varid,varidx,varidy,len_file,err
  real, dimension(:,:,:,:), allocatable :: swd,A,T
  real, dimension(:,:), allocatable :: x,y
  integer, allocatable, dimension(:,:,:) :: dims
  integer :: xdim,ydim,tdim,NY

  ! integer,parameter :: xdim=192,ydim=96,tdim=12,NY=1
  ! real  :: x(xdim,1)
  ! real  :: y(ydim,1)
  ! ! real  :: lon(LLm,MMm),lat(LLm,MMm)
  ! real  :: swd(xdim,ydim,tdim,NY),A(xdim,ydim,tdim,NY),T(xdim,ydim,tdim,NY)
  !
  character*300 :: file_in

  file_in='/Users/xushan/Downloads/dEBM-master/pi477_u_MM_2601_270001.01_echam.nc'
  len_file=len_trim(file_in)
  ierr=nf_open(trim(file_in),nf_nowrite,ncid)
  ierr=nf_inq_varid(ncid,'rsuscs_na',varid)
  ierr=nf_get_var_real(ncid,varid,swd)

  ! dims = shape(swd)
  write(*,*) shape(swd)

 write(*,*) swd(159,45,1,1),swd(154,4,1,1)


  ! xdim   = size(swd[:], 1)
  ! ydim   = size(swd[:], 2)
  ! tdim   = size(swd[:], 3)
  ! NY     = size(swd[:], 4)

  write(*,*) xdim,ydim,tdim,NY

  !allocate (x(xdim,ydim))
  !allocate (y(xdim,ydim))
  !allocate (T(xdim,ydim,tdim,NY))
  !allocate (A(xdim,ydim,tdim,NY))
  !allocate (swd(xdim,ydim,tdim,NY))
  !
  ! ierr=nf_inq_varid(ncid,'temp2',varid)
  ! ierr=nf_get_var_real(ncid,varid,T)
  ! ierr=nf_inq_varid(ncid,'alsoi',varid)
  ! ierr=nf_get_var_real(ncid,varid,A)
  ! ierr=nf_inq_varid(ncid,'lon',varidx)
  ! ierr=nf_inq_varid(ncid,'lat',varidy)
  ! ierr=nf_get_var_real(ncid,varidx,x)
  ! ierr=nf_get_var_real(ncid,varidy,y)
  !
  !print *,swd(159,45,1,1),swd(154,4,1,1)
  err=NF_CLOSE(ncid)
END PROGRAM TEST
