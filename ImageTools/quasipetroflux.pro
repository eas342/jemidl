Function QuasiPetroFlux, values, eta=eta
  if n_elements(eta) eq 0 then eta=0.2
  s=reverse(sort(values))
  cumulative = total( values[s], /cumulative )
  for i=0L, n_elements(values)-1 do begin
     if values[s[i]] - eta * cumulative[i]/(i+1) lt 0 then break
  endfor
  if i eq n_elements(values) then return, !values.f_nan
  return, values[s[i]]
end
