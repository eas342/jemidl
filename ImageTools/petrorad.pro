Function PetroRad, image, center, boa, pa, eta=eta
  if n_elements(eta) eq 0 then eta=0.2
  s=size(image,/dim)
  maxa = max(s)/2
  a = findgen(maxa)+1
  flux = apellipphot( image, center[0]+0.5, center[1]+0.5, boa, pa, a )
  fluxdens = flux/(a*a)
  a2 = [0,a[0:n_elements(a)-2]]
  flux2 = [0,flux[0:n_elements(flux)-2]]
  anfluxdens = ((flux-flux2)/((a^2-a2^2))[0:n_elements(a)-2])
  y = eta*fluxdens-anfluxdens
  if y[0]*max(y) gt 0. then return, !values.f_nan
  return, (interpol(0.5*(a[0:n_elements(a)-2]+a[1:*]), y, [0.]))[0]
end
