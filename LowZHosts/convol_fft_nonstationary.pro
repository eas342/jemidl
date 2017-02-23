Function convol_fft_nonstationary, f, kernel
  d=size(f, /n_dim)
  if d eq 1 then begin
     nfunc=1
     npix=(size(f,/dim))[0]
  endif else begin
     s=size(f, /dim)
     npix = s[0]
     nfunc = s[1]
  endelse
  d=size(kernel, /n_dim)
  ;; if d eq 1 then begin ;; specified as sigma
  ;;    npixel=2*ceil(4*max(kernel))+1
  ;;    kernel1 = dblarr(npixel, npix)
  ;;    for i=0, npix-1 do begin
  ;;       kernel1[*,i] = psf_Gaussian(npixel=npixel, st_dev=kernel[i], /norm, ndim=1)
  ;;    endfor
  ;; endif else begin
  ;;    kernel1 = kernel
  ;; endelse
  out = f*0
  for i=0L, npix-1 do begin
     if d eq 1 then begin
        npixel=2*ceil(4*kernel[i])+1
        kernel1 = psf_Gaussian(npixel=npixel, st_dev=kernel[i], /norm, ndim=1)
     endif else begin
        kernel1 = kernel[*,i]
        npixel=n_elements(kernel1)
     endelse
     ;; to speed things up, only convolve region near pixel of interest
     lo = i-npixel-1
     hi = i+npixel+1
     if lo lt 0 then begin
        lo=0
        hi=2*npixel+2
     endif
     if hi gt npix-1 then begin
        hi=npix-1
        lo=npix-2*npixel-3
     endif
     out[i,*] = (convol_fft2(f[lo:hi,*], kernel1))[i-lo,*]
  endfor
  return, out
end
