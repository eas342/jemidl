Function Sersic_flight, n, r
  k = Sersic_Kappa(n)
  return, igamma(2*n, r^(1/n)*k)
end
