Function QuasiPetroRad, image, mask, center, boa, pa, eta=eta, qpflux=qpflux
  qpflux = quasipetroflux(image*mask)
  if n_elements(eta) eq 0 then eta=0.2
  s=size(image,/dim)
  maxa = max(s)/2
  a = findgen(maxa)+1
  flux = apellipphot( image, center[0]+0.5, center[1]+0.5, boa, pa, a)
  fluxdens = flux/(!dpi*a*a*boa)*eta
  return, (interpol(a, fluxdens, [qpflux]))[0]
end
