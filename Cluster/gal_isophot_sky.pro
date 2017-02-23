Pro Gal_isophot_sky, img_ptr, err_ptr, gal, galapagoscat, platescale=platescale, sky_value=sky_value, sky_sigma=sky_sigma, sky_median=sky_median, skyim=skyim, skyreg=skyreg

  if n_elements(platescale) eq 0 then platescale = 0.05
  platefrac = 0.05/platescale
  ;; first get an extended postage stamp...
  xmin = gal.xmin - 250*platefrac > 0
  xmax = gal.xmax + 250*platefrac < (size(*img_ptr, /dim))[0]-1
  ymin = gal.ymin - 250*platefrac > 0
  ymax = gal.ymax + 250*platefrac < (size(*img_ptr, /dim))[1]-1

  nx = xmax - xmin + 1
  ny = ymax - ymin + 1

  skyim = (*img_ptr)[xmin:xmax, ymin:ymax]
  skyer = (*err_ptr)[xmin:xmax, ymin:ymax]
  mask = skyim*0
  
  ;;center is off the edge...
  if ~finite((*err_ptr)[(gal.xmin+gal.xmax)/2,(gal.ymin+gal.ymax)/2]) then begin
     sky_value=!values.d_nan
     sky_sigma=!values.d_nan
     sky_median=!values.d_nan
     return
  endif
  ;; make ellipses
  ellipses=bytarr( nx, ny, 14 )
  annuli  =bytarr( nx, ny, 12 )
  dist_ellipse, skyellipse, [nx,ny], $
                gal.x_guess-xmin-1, gal.y_guess-ymin-1, $
                1./gal.boa_guess, gal.pa_guess
  for iellipse=0, 13 do begin
     ellipses[*,*,iellipse] = skyellipse lt gal.majaxis + 18*iellipse*platefrac
  endfor
  
  for iannulus=0, 11 do begin
     annuli[*,*,iannulus] = ellipses[*,*,iannulus+2] and (1 - ellipses[*,*,iannulus])
  endfor
  

  ;; mask out secondary objects
  w = where( galapagoscat.xmin le xmax and $
             galapagoscat.xmax ge xmin and $
             galapagoscat.ymin le ymax and $
             galapagoscat.ymax ge ymin )
  if w[0] ne -1 then begin
     for j=0, n_elements(w)-1 do begin
        maskgal = galapagoscat[w[j]]
        if maskgal.galid eq gal.galid then continue
        dist_ellipse, maskellipse, [nx,ny], $
                      maskgal.x_guess-xmin-1, maskgal.y_guess-ymin-1, $
                      1./maskgal.boa_guess, maskgal.pa_guess
        e = maskellipse le maskgal.majaxis
        mask = mask or e
     endfor
  endif
  
  ;; mask out unexposed pixels
  mask = mask or (1-finite( skyer ))
  
  ;; calculate mean sky value in each annulus
  sky = dblarr(12)
  minstart=0
  for iannulus=0, 11 do begin
     wsky = where( annuli[*,*,iannulus] and (1-mask) )
     if n_elements(wsky) lt 20 then minstart += 1 $
     else sky[iannulus] = biweight_mean( skyim[wsky] )
  endfor
  
  grad=dblarr(7)
  for igrad=0,6 do begin
     result = linfit(indgen(6), sky[igrad:igrad+5])
     grad[igrad]=result[1]
  endfor

  wgrad = where(grad gt -0.003/(18*platefrac)*0.01)
  if wgrad[0] eq -1 then istart=6 $
  else istart=max([min(wgrad),minstart]) < 6
;  atv, skyim*(ellipses[*,*,istart+6] and (1-ellipses[*,*,istart]) and (1-mask))
  w = where( ellipses[*,*,istart+6] and (1-ellipses[*,*,istart]) and (1-mask) )
  if w[0] eq -1 then begin
     sky_value=!values.d_nan
     sky_median=!values.d_nan
     sky_sigma=!values.d_nan
     return
  endif
  if arg_present(skyreg) then skyreg=skyim*(ellipses[*,*,istart+6] and (1-ellipses[*,*,istart]) and (1-mask))
  sky_value = biweight_mean(skyim[w], sky_sigma )
  sky_median = median(skyim[w])
  return
end
