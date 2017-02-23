Function cheb_eval, coeff, x, domain=domain
  if n_elements(domain) eq 0 then domain=minmax(x)
  n=n_elements(coeff)
  array = cheb_array(n-1, x, domain=domain)
  return, array#coeff
end
