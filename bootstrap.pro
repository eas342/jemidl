Function bootstrap, array
  m = fltarr(1000)
  s = fltarr(1000)
  for i=0,999 do begin
     w=floor(randomu(seed, n_elements(array))*n_elements(array))
     m[i] = biweight_mean(array[w], sigma)
     s[i] = sigma
  endfor
  return, [mean(m), stdev(m), mean(s), stdev(s)]
end
