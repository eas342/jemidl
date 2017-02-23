Function SchechterIntegrand, M
  common JEM$_SchechterIntegrand, phi_star1, M_star1, alpha1
  return, SchechterM( M, phi_star1, M_star1, alpha1 )
end

Function SchechterMBin, phi_star, M_star, alpha, Mmin, Mmax
  common JEM$_SchechterIntegrand
  phi_star1=phi_star
  M_star1=M_star
  alpha1=alpha
  return, qpint1d( 'SchechterIntegrand', Mmin, Mmax )
end
