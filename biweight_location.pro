Function MedAbsDev, $
   x, $
   M

  if n_elements(M) eq 0 then M=median(x)
  return, median(abs(x-M))
end

Function biweight_location1, $
   x, $
   M, $
   MAD

  u = (x - M) / (6.0 * MAD)
  w = where(abs(u) lt 1.0)
  if w[0] eq -1 then return, M
  CBI = M + total( (x - M) * (1-u[w]^2)^2) / total((1-u[w]^2)^2)
  return, CBI
end

Function biweight_location, $
   x, $
   nan=nan

  if keyword_set(nan) then begin
     w=where(finite(x))
     if w[0] eq -1 then return, !values.f_nan $
     else x = x[w]
  endif
  CBI = median(x)
  for i=0, 3 do begin
     M = CBI
     MAD = MedAbsDev( x, M )
     CBI = biweight_location1( x, M, MAD )
  endfor
  return, CBI
end
