;josh sky aperture errors
Function geterrs5, obj
  id=obj->extract('clusterid')
  restore, '/home/scpdata02/joshimages/skystats/skysummary.sav'
  w=where(skystats.id eq id)
  s=obj->summary()
  c=obj->color(/xpsf, mags=mags)
  zp=obj->extract('zeropoint')
  r=s.re > 3.
  ap = [1.0, 1.3, 1.7, 2.2, 2.8, 3.6, 4.6, 6.0, 7.7, 10.0]
  iline=linfit( ap, skystats[w].icvlsig)
  zline=linfit( ap, skystats[w].zcvlsig)

  wsmall=where(r le 10., complement=wbig)
  ifluxerr=fltarr(n_elements(c))
  zfluxerr=fltarr(n_elements(c))

  ifluxerr[wsmall] = interpol(skystats[w].icvlsig, ap, r[wsmall])
  zfluxerr[wsmall] = interpol(skystats[w].zcvlsig, ap, r[wsmall])
  ifluxerr[wbig] = iline[0]+iline[1]*r[wbig]
  zfluxerr[wbig] = zline[0]+zline[1]*r[wbig]
  iflux = 10^(-0.4*(mags[*,0]-zp[0]))
  zflux = 10^(-0.4*(mags[*,1]-zp[1]))
  iexptime = obj->sxpar('EXPTIME', band='i')
  zexptime = obj->sxpar('EXPTIME', band='z')
  ifluxerr = sqrt(ifluxerr^2+iflux/iexptime)
  zfluxerr = sqrt(zfluxerr^2+zflux/zexptime)
  idfof = ifluxerr/iflux
  zdfof = zfluxerr/zflux
  ierr = 2.5/2*alog10((1+idfof)/(1-idfof))
  zerr = 2.5/2*alog10((1+zdfof)/(1-zdfof))
  return, sqrt(ierr^2+zerr^2)
end
