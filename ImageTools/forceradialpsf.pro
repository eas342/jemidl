Function ForceRadialPSF, psf, set, bkspace=bkspace, loc=loc, scale=scale
  origset = GetRadialPSF( psf, bkspace=bkspace )
  s=size(psf,/dim)
  y=alog10(abs(psf)/max(psf, /nan))
  dist_ellipse, dist, s, (s[0]-1)/2, (s[1]-1)/2, 1, 0
  w=where(finite(y))
  targ = bspline_valu( dist[w], set )
  orig = bspline_valu( dist[w], origset )
  activate = 1./!pi*atan(scale*(dist[w]-loc))+0.5
;  y[w] *= 1.+(targ-orig)*activate/orig
  y[w] += (targ-orig)*activate
  return, 10.^y
end
