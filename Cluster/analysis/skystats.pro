Function rmelt, array, elts
  out = array
  for ielt=0, n_elements(elts)-1 do begin
     w=where(out eq elts[ielt], complement=wkeep)
     if wkeep[0] eq -1 then return, -1
     out = out[wkeep]
  endfor
  return, out
end

Function skystats, image, cvlimage, mask, apertures, niter=niter
  if n_elements(niter) eq 0 then niter=1000
  nap = n_elements(apertures)
  maxap=max(apertures)
  tmpmask = mask
  for ishift=-maxap, maxap do begin
     for jshift=-maxap, maxap do begin
        counter, (2*maxap+1)*(ishift+maxap)+(jshift+maxap+1), (2*maxap+1)^2
        if ishift^2+jshift^2 gt maxap^2 then continue
        tmpmask = tmpmask and shift(mask, ishift, jshift)
     endfor
  endfor

  oneskystat = { center: fltarr(2), $
                 fluxes: fltarr(nap), $
                 cvlfluxes: fltarr(nap), $
                 sky_value: 0., $
                 sky_sigma: 0. $
               }
  skystats = replicate(oneskystat, niter)
  sz = size(image, /dim)
  for i=0, niter-1 do begin
     counter, i+1, niter
     wgood=where(tmpmask)
     if wgood[0] eq -1 then break
     p=randomu(seed)*n_elements(wgood)
     center=array_indices(tmpmask, wgood[p])
     center += randomu(seed, 2)
     skyxmin = floor(center[0]-250 > 0)
     skyymin = floor(center[1]-250 > 0)
     skyxmax = floor(center[0]+250 < sz[0]-1)
     skyymax = floor(center[1]+250 < sz[1]-1)
     ellipse = { ellipse, $
                 center: center-[skyxmin,skyymin], $
                 majaxis: maxap, $
                 minaxis: maxap, $
                 theta: 0.}
     skystruct = newisophotsky( image[skyxmin:skyxmax, skyymin:skyymax], $
                                ~tmpmask[skyxmin:skyxmax, skyymin:skyymax], $
                                ellipse, $
                                platescale=0.05, skyreg=skyreg )
     xmin = floor(center[0] - 1.5*max(apertures) > 0 )
     xmax = floor(center[0] + 1.5*max(apertures) < sz[0]-1 )
     ymin = floor(center[1] - 1.5*max(apertures) > 0 )
     ymax = floor(center[1] + 1.5*max(apertures) < sz[1]-1 )
     
     dist_ellipse, minimask, [xmax-xmin+1, ymax-ymin+1], $
                   center[0]-xmin, center[1]-ymin, 1., 0.
     minimask = minimask ge 2*max(apertures)+3
     tmpmask[xmin:xmax, ymin:ymax] = $
        tmpmask[xmin:xmax, ymin:ymax] and minimask

     subimage = image[xmin:xmax, ymin:ymax] - skystruct.sky_value
     subcvlimage = cvlimage[xmin:xmax, ymin:ymax] - skystruct.sky_value
     flux1 = apphot( subimage, $
                     center[0]-xmin+randomu(seed), $
                     center[1]-ymin+randomu(seed), $
                     apertures )
     flux2 = apphot( subcvlimage, $
                     center[0]-xmin+randomu(seed), $
                     center[1]-ymin+randomu(seed), $
                     apertures )
     skystats[i].fluxes=flux1
     skystats[i].cvlfluxes=flux2
     skystats[i].center=center
     skystats[i].sky_value=skystruct.sky_value
     skystats[i].sky_sigma=skystruct.sky_sigma
  endfor
  return, skystats[0:i-1]
end
