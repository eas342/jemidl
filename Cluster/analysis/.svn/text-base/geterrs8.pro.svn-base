Function geterrs8, $
   re, $
   isky_sigma, $
   zsky_sigma, $
   imag, $
   zmag, $
   iexptime, $
   zexptime, $
   ierr=ierr, $
   zerr=zerr, $
   zp=zp

  ngal = n_elements(re)
  re = (finite(s.re)*s.re+1.2 > 3)*finite(s.re) ;; take max( re+1.2, 3.0 ) as re
  wsmall=where(re le 10.)
  wbig=where(re gt 10.)

  ifluxerrs = dblarr(ngal)
  zfluxerrs = dblarr(ngal)

  if wsmall[0] ne -1 then begin
     for iobj=0, n_elements(wsmall)-1 do begin
        ifluxerrs[wsmall[iobj]] = $
           interpol( ifits[*,0]+ifits[*,1]*s[wsmall[iobj]].isky_sigma, $
                     radii, $
                     re[wsmall[iobj]] )
        zfluxerrs[wsmall[iobj]] = $
           interpol( zfits[*,0]+zfits[*,1]*s[wsmall[iobj]].zsky_sigma, $
                     radii, $
                     re[wsmall[iobj]] )
     endfor
  endif

  if wbig[0] ne -1 then begin
     for iobj=0, n_elements(wbig)-1 do begin
        fit = linfit( radii[7:*], ifits[7:*,0]+ifits[7:*,1]*s[wbig[iobj]].isky_sigma )
        ifluxerrs[wbig[iobj]] = fit[0]+fit[1]*re[wbig[iobj]]
        fit = linfit( radii[7:*], zfits[7:*,0]+zfits[7:*,1]*s[wbig[iobj]].zsky_sigma )
        zfluxerrs[wbig[iobj]] = fit[0]+fit[1]*re[wbig[iobj]]
     endfor
  endif

  if n_elements(zp) eq 0 then zp = [24.8666, 25.6785]
  ifluxes = 10^(-0.4*(imag-zp[1]))
  zfluxes = 10^(-0.4*(zmag-zp[0]))
  ifluxerrs = sqrt(ifluxerrs^2 + ifluxes/iexptime)
  zfluxerrs = sqrt(zfluxerrs^2 + zfluxes/zexptime)
  idfof = ifluxerrs/ifluxes
  zdfof = zfluxerrs/zfluxes
  ierr = 2.5/2*alog10((1+idfof)/(1-idfof))
  zerr = 2.5/2*alog10((1+zdfof)/(1-zdfof))
  return, sqrt(ierr^2+zerr^2)
end
