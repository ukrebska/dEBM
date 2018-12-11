! again just for optimization
subroutine dEBM_cosphicosd(lat,decl)
real :: lat,decl
real, parameter :: pi=3.1425926
cosphicosd=COS(pi/180*lat).*COS(pi/180*decl)
end subroutine dEBM_cosphicosd
