Function powerlist, list, factor, max
  powermax = floor( alog(max) / alog(factor) )
  powers = ulon64arr(powermax+1)
  powers[0] = 1
  for i=1, powermax do powers[i] = factor*powers[i-1]
  listmat = list#powers
  inds = indgen(n_elements(powers))
  maxind = floor( (alog(max/list)) / alog(factor) )
  for i=0, n_elements(list)-1 do begin
     w=where(inds gt maxind[i],count)
     if count gt 0 then listmat[i,w] = 1
  endfor
  return, listmat[uniq(listmat, sort(listmat))]
end

Function powers, factors, max
  list = [1ULL]
  for i=0, n_elements(factors)-1 do begin
     list = powerlist( list, factors[i], max )
  endfor
  return, list
end
