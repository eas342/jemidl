;; stolen from ppxf
Function convol_fft2, f, k
  d=size(f, /n_dim)
  if d eq 1 then begin
     nfunc=1
     npix=(size(f,/dim))[0]
  endif else begin
     s=size(f, /dim)
     npix = s[0]
     nfunc = s[1]
  endelse
  out = f

  nf = npix
  nk = n_elements(k)
  n = 2L^ceil(alog(nf+nk/2)/alog(2))
  for i=0, nfunc-1 do begin
     f1 = dblarr(n)
     k1 = f1
     f1[0] = f[*,i]
     k1[0] = rotate(k,2)
     k1 = shift(k1,-(nk-1)/2)
     con = n*double(fft(fft(f1,-1)*fft(k1,-1),1))
     out[*,i]=con[0:nf-1]
  endfor
  return, out
end
