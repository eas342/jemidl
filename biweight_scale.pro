Function MedAbsDev, $
   x, $
   M

  if n_elements(M) eq 0 then M=median(x)
  return, median(abs(x-M))
end

Function biweight_scale1, $
   x, $
   M, $
   MAD, $
   u=u, $
   c=c

  u = (x - M) / (c * MAD)
  w = where(abs(u) lt 1.0)
  if w[0] eq -1 then return, !values.f_infinity
  SBI = sqrt(n_elements(x))*sqrt(total((x[w]-M)^2*(1-u[w]^2)^4)) $
             / abs(total((1-u[w]^2)*(1-5*u[w]^2)))
  return, SBI
end

Function biweight_scale, $
   x, $
   zero=zero, $
   weight=weight, $
   c=c, $
   MAD=MAD, $
   niter=niter, $
   nan=nan

  if keyword_set(nan) then begin
     w=where(finite(x))
     if w[0] eq -1 then return, !values.f_nan $
     else x = x[w]
  endif
  if n_elements(niter) eq 0 then niter=4
  if n_elements(c) eq 0 then c=9.0
  if keyword_set(zero) then M = 0 $
  else M = biweight_location(x)
  if n_elements(MAD) eq 0 then MAD = MedAbsDev( x, M )
  SBI = MAD/0.6745
  for i=0, niter-1 do begin
     MAD = 0.6745*SBI
     SBI = biweight_scale1( x, M, MAD, u=u, c=c )
  endfor
  weight=1-(u^2.<1.)
  return, SBI
end
