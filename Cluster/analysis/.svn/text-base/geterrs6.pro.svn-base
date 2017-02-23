;sky sigma errorbars
Function geterrs6, obj
  xpsfap = [3.0, 5.0, 10.0, 15.0, 20.0]
  id=obj->extract('clusterid')
  s=obj->summary()
  c=obj->color(/xpsf, mags=mags)
  re = s.re > 3
  A = sqrt(!pi)*re
  ifluxerrs = s.isky_sigma*A*9./7
  zfluxerrs = s.zsky_sigma*A*9./7
  zp=obj->extract('zeropoint')
  ifluxes = 10^(-0.4*(mags[*,1]-zp[1]))
  zfluxes = 10^(-0.4*(mags[*,0]-zp[0]))
  idfof = ifluxerrs/ifluxes
  zdfof = zfluxerrs/zfluxes
  ierr = 2.5/2*alog10((1+idfof)/(1-idfof))
  zerr = 2.5/2*alog10((1+zdfof)/(1-zdfof))
  return, sqrt(ierr^2+zerr^2)
end
