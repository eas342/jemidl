Function sampleskewness, X
  n=n_elements(X)*1LL
  mu = total(X)/n
  sig2 = total(X^2)/n - mu^2
  m3 = (total(X^3)/n - 3*mu*sig2 - mu^3)
  g1 = m3/(sig2^1.5)
  return, sqrt(n*(n-1))/(n-2)*g1
end
