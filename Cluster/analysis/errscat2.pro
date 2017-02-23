Function errscat2, errs
  scat = fltarr(5000)
  for i=0, 4999 do begin
     arr = randomn(seed, n_elements(errs))*errs
;     scat[i] = robust_sigma(arr, /zero)
     scat[i] = biweight_scale( arr, /zero )
  endfor
  return, mean(scat)
end
