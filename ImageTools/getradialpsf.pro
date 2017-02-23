Function GetRadialPSF, psf, bkspace=bkspace
  s=size(psf,/dim)
  y=alog10(abs(psf)/max(psf, /nan))
  dist_ellipse, dist, s, (s[0]-1)/2, (s[1]-1)/2, 1, 0
  s=sort(dist)
  x=dist[s]
  y=y[s]
  w=where(finite(y))
  x=x[w]
  y=y[w]
  fullbkpt = bspline_bkpts( x, bkspace=bkspace, nord=4 )
  err = bspline_fit( x, y, fltarr(n_elements(y))+1., set, fullbkpt=fullbkpt )
  return, set
end
