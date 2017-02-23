;first principles errorbars with correction for correlated pixels
Function geterrs4, obj
  xpsfap = [3.0, 5.0, 10.0, 15.0, 20.0]
  id=obj->extract('clusterid')
  s=obj->summary()
  c=obj->color(/xpsf, mags=mags)
  xpsfcat = obj->extract('xpsfcat')
  re = s.re > 3
  tabinv, xpsfap, re, apindex
  ierr=fltarr(n_elements(s))
  zerr=fltarr(n_elements(s))
  for igal=0, n_elements(s)-1 do begin
     ierr[igal] = interpolate( transpose(xpsfcat[igal,1].xpsfaperr), apindex[igal] )
     zerr[igal] = interpolate( transpose(xpsfcat[igal,0].xpsfaperr), apindex[igal] )
  endfor
  return, sqrt(ierr^2+zerr^2)
end
