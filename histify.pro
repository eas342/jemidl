Function histify, x, y, fill=fill
  bin = x[1]-x[0]
  nbin = n_elements(x)
  binleft = x - bin/2.0
  binright = x + bin/2.0
  xout = fltarr(nbin*2+2+keyword_set(fill))
  yout = fltarr(nbin*2+2+keyword_set(fill))
  xout[0] = binleft[0]
  yout[0] = 0.0
  for i=0, nbin-1 do begin
     xout[2*i+1] = binleft[i]
     xout[2*i+2] = binright[i]
     yout[2*i+1] = y[i]
     yout[2*i+2] = y[i]
  endfor
  xout[2*nbin+1] = binright[nbin-1]
  yout[2*nbin+1] = 0.0

  if keyword_set(fill) then begin
     xout[2*nbin+2]=binleft[0]
     yout[2*nbin+2]=0.0
  endif

  return, [[xout], [yout]]
end
