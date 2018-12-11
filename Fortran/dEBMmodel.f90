subroutine dEBMmodel(A, swd, PDD,hours,q,c1,c2,MC)
implicit none

real, dimension(:,:,:,:) :: A,swd,PDD,hous,q,MC,melt
real, paramter :: Lf=3.34e5         ! Lf is latent heat of fusion

! melt is the melt rate (mm/day)
dEBMmelt=MC.*max(0,((1-A).*swd.*24.*q)+(c2+c1*PDD).*hours)/Lf*60*60
return

end subroutine dEBMmodel
