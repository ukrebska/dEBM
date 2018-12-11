! if F is the flux density normal to a surface,
! insolation between [-ha,ha] is SW=  sol_flux_fact*F
! ha is an hour angle during day time
subroutine dEBM_sol_flux_fact(sin_phi_sin_d,cos_phi_cos_d,ha)

solfluxfact=(sin_phi_sin_d.*ha+cos_phi_cos_d.*sin(ha))/pi

end subroutine dEBM_sol_flux_fact
