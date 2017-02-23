;skysig => skyap errors
Function geterrs7, obj, ierr=ierr, zerr=zerr
  restore, '/home/scpdata02/joshimages/skystats/sky2apsig.sav'
  s=obj->summary()
  junk=obj->color(/xpsf, mags=mags, /silent)
  re = (finite(s.re)*s.re+1.2 > 3)*finite(s.re)
  wsmall=where(re le 10.)
  wbig=where(re gt 10.)
  ifluxerrs = dblarr(n_elements(s))
  zfluxerrs = dblarr(n_elements(s))
  for iobj=0, n_elements(wsmall)-1 do begin
     ifluxerrs[wsmall[iobj]] = interpol( ifits[*,0]+ifits[*,1]*s[wsmall[iobj]].isky_sigma, radii, re[wsmall[iobj]] )
     zfluxerrs[wsmall[iobj]] = interpol( zfits[*,0]+zfits[*,1]*s[wsmall[iobj]].zsky_sigma, radii, re[wsmall[iobj]] )
  endfor
  for iobj=0, n_elements(wbig)-1 do begin
     fit = linfit( radii[7:*], ifits[7:*,0]+ifits[7:*,1]*s[wbig[iobj]].isky_sigma )
     ifluxerrs[wbig[iobj]] = fit[0]+fit[1]*re[wbig[iobj]]
     fit = linfit( radii[7:*], zfits[7:*,0]+zfits[7:*,1]*s[wbig[iobj]].zsky_sigma )
     zfluxerrs[wbig[iobj]] = fit[0]+fit[1]*re[wbig[iobj]]
  endfor
  zp=obj->extract('zeropoint')
  ifluxes = 10^(-0.4*(mags[*,1]-zp[1]))
  zfluxes = 10^(-0.4*(mags[*,0]-zp[0]))
  ifluxerrs = sqrt(ifluxerrs^2 + ifluxes/s.iexptime)
  zfluxerrs = sqrt(zfluxerrs^2 + zfluxes/s.zexptime)
  idfof = ifluxerrs/ifluxes
  zdfof = zfluxerrs/zfluxes
  ierr = 2.5/2*alog10((1+idfof)/(1-idfof))
  zerr = 2.5/2*alog10((1+zdfof)/(1-zdfof))
  return, sqrt(ierr^2+zerr^2)
end
