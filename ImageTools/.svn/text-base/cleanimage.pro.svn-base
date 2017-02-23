Function CLEANImage, image, psf, fwhm, outims=outims, maxiter=maxiter
  if n_elements(maxiter) eq 0 then maxiter=500
  xs = indgen(7)+0.5
  ys = indgen(7)+0.5
  sig = fwhm/(2*sqrt(2*alog(2)))
  param = [0,1,sig,sig,3.5,3.5]
  psf /= total(psf)
  gauss = sampledgauss2d( xs, ys, param )
  gauss /= total(gauss)

  if arg_present(outims) then begin
     s=size(image,/dim)
     outims=fltarr(s[0]*3,s[1],maxiter)
  endif

  if max(image) le 0. then return, image
;  cln = clean3( image, psf, 0.001, 2000 )
  cln = clean3( image, psf, 0.001, 200 )
  gain = 0.01
  for i=0, maxiter-1 do begin
     nextcln = clean3( cln[*,*,1], psf, gain, 2000 )
     cln[*,*,0] += nextcln[*,*,0]
     cln[*,*,1]  = nextcln[*,*,1]

     max = max(cln[*,*,1])
     junk=biweight_mean(cln[*,*,1], sigma)
     s2n = max/sigma

     if s2n gt 15. then gain=0.01 $
     else if s2n gt 7. then gain=0.003 $
     else if s2n gt 5. then gain=0.001 $
     else break
     
     if arg_present(outims) then begin
        outims[*,*,i]=[convol(cln[*,*,0],gauss,/edge_zero), $
                       cln[*,*,1], $
                       convol(cln[*,*,0],gauss,/edge_zero)+cln[*,*,1]]
     endif
  endfor
  return, convol(cln[*,*,0],gauss,/edge_zero)+cln[*,*,1]
end
