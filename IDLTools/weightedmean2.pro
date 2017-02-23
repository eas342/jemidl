Function weightedmean2test, sigint1
  common weightedmean2block, x, w, m, sigint
  result = fltarr(n_elements(sigint1))
  for i=0, n_elements(result)-1 do begin
     result[i] = total((x-m)^2/(1/w+sigint1[i]^2))-(n_elements(x)-2)
  endfor
  return, result
end

Function weightedmean2test2, mean
  common weightedmean2block
  result = fltarr(n_elements(mean))
  for i=0, n_elements(result)-1 do begin
     result[i] = total((x-mean[i])^2/(1/w+sigint^2))-(n_elements(x)-2)
  endfor
  return, result
end

Function weightedmean2, x1, w1, sigint=sigint1
  common weightedmean2block
  x = x1
  w = w1
  m = weightedmean(x, w)
  rms = stdev(x)
  sigsqr = mean(1/w1)
  sigintestimate = sqrt(rms^2-sigsqr)
  chi2 = total((x-m)^2*w)
  if chi2 lt n_elements(x)-2 then begin
     xsigma=!values.f_nan
     return, m
  endif
  sigint1 = fx_root([sigintestimate*0.9, sigintestimate, sigintestimate*1.1], 'weightedmean2test')
  sigint=sigint1
  stop
;  sigint=sigint1
;  up = fx_root([m,m*1.01,m*1.02], 'weightedmean2test2')
;  down = fx_root([m*0.98,m*0.99,m], 'weightedmean2test2')
  return, m
end
