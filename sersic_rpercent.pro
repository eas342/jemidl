Function igammakappa, kappa
  common sersic_block, n1, f1
  return, igamma(2*n1, kappa)-f1
end

Function Sersic_RPercent, n, percent
  common sersic_block
  k = Sersic_Kappa(n)
  f1 = percent/100.0
  alpha = newton(k, 'igammakappa', stepmax=0.1)
  return, (alpha/k)^n
end
