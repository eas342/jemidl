Function errscat, errs
  scat=fltarr(5000)
  for i=0, 4999 do begin
     arr = randomn(seed, n_elements(errs))*errs
     junk=biweight_mean(arr, sigma)
     scat[i]=sigma
  endfor
  return, mean(scat)
end
