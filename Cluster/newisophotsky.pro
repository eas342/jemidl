Function newisophotsky, image, mask, objellipse, platescale=platescale, skyreg=skyreg
  if n_elements(platescale) eq 0 then platescale=0.05
  platefrac = 0.05/platescale
  sz = size( image, /dim )
  ellipses = bytarr( sz[0], sz[1], 14 )
  annuli   = bytarr( sz[0], sz[1], 12 )

  @definestructs
  out=isophotsky0

  dist_ellipse, skyellipse, sz, $
                objellipse.center[0], objellipse.center[1], $
                objellipse.majaxis/objellipse.minaxis, objellipse.theta-90.
  for iellipse=0, 13 do begin
     ellipses[*,*,iellipse] = skyellipse lt objellipse.majaxis+18*iellipse*platefrac
  endfor
  for iannulus=0, 11 do begin
     annuli[*,*,iannulus] = ellipses[*,*,iannulus+2] and (1 - ellipses[*,*,iannulus])
  endfor
  
  sky = dblarr(12)
  minstart=0
  for iannulus=0, 11 do begin
     wsky = where( annuli[*,*,iannulus] and (1-mask) )
     if n_elements(wsky) lt 20 then minstart += 1 $
     else sky[iannulus] = biweight_mean( image[wsky] )
  endfor

  grad=dblarr(7)
  for igrad=0,6 do begin
     result = linfit(indgen(6), sky[igrad:igrad+5])
     grad[igrad]=result[1]
  endfor
  wgrad = where(grad gt -0.003/(18*platefrac)*0.01)
  if wgrad[0] eq -1 then istart=6 $
  else istart=max([min(wgrad), minstart]) < 6
  w = where( ellipses[*,*,istart+6] and (1-ellipses[*,*,istart]) and (1-mask) )
  if w[0] ne -1 then begin
     if arg_present(skyreg) then skyreg=image*(ellipses[*,*,istart+6] $
                                               and (1-ellipses[*,*,istart]) $
                                               and (1-mask))
     out.sky_value = biweight_mean( image[w], sigma )
     out.sky_sigma = sigma
     out.sky_median = median(image[w])
  endif
  return, out
end
