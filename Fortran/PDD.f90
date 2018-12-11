subroutine PDD(T,stddev)
! PDD is approximated as described in
! Calov R and Greve R (2005)
! Correspondence. A semi-analytical solution for the positive degree-day model
! with stochastic temperature variations.
! J. Glaciol., 51 (172), 173–175
! (doi:10.3189/172756505781829601)
!****************************************************************************
! T in ^oC!!!
! stddev can be spatially variable (same grid as T) or a constant.
!T=T-273.15;

PDD=stddev/sqrt(2*pi)*exp(-(T**2)/(2*(stddev**2)))+T/2.*erfc(-T/(sqrt(2)*stddev))!
end subroutine PDD
