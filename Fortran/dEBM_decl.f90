subroutine dEBM_decl(days,obliqu)
decl=obliqu*sin(pi/180*360/365*(days(:)-79));   ! 79 is the number of days since New Year for the spring equinox
end subroutine dEBM_decl
