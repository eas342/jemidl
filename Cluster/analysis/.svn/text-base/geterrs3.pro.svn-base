;kyle barbary's skystats.
Function geterrs3, obj, _extra=extra
  id=obj->extract('clusterid')
  restore, '/home/scpdata02/clusters/skystats.sav'
  s=obj->summary()
  c=obj->color(/xpsf, mags=mags)
  zp=obj->extract('zeropoint')
  r=s.re > 3.

  if id eq 'D' then id='O'
  wi=where(where(strpos(filename, 'CL-'+id) ne -1) and $
           strpos(filename, 'F775W') ne -1)
  wz=where(where(strpos(filename, 'CL-'+id) ne -1) and $
           strpos(filename, 'F850LP') ne -1)

  iline=linfit( ap[wi], sky[wi] )
  zline=linfit( ap[wz], sky[wz] )
  wsmall=where(r le 10., complement=wbig)
  ifluxerr=fltarr(n_elements(c))
  zfluxerr=fltarr(n_elements(c))

  ifluxerr[wsmall] = interpol(sky[wi], ap[wi], r[wsmall])
  zfluxerr[wsmall] = interpol(sky[wz], ap[wz], r[wsmall])
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
