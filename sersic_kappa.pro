Function igammakappa, kappa
  common sersic_block, n1, f1
  return, igamma(2*n1, kappa)-f1
end

Function Sersic_Kappa, n
  common sersic_block
  n1 = n
  f1 = 0.5
  kappaguess = 2*n - (1./3)
  return, newton(kappaguess, 'igammakappa', stepmax = 0.1)
end

