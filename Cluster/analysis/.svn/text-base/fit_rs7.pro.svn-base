Pro fit_rs7, $
   cmd, $
   slope, $
   magrange=magrange, $
   bluelimit=bluelimit, $
   diagnostic=diagnostic

  if n_elements(bluelimit) eq 0 then bluelimit=0.
  if n_elements(magrange) eq 0 then magrange=[19.,25.]
  gals=*cmd.gals
  wselect = where( gals.select $
                   and gals.z850 ge magrange[0] $
                   and gals.z850 le magrange[1] $
                   and gals.iz ge bluelimit )
  select = gals[wselect]
  z850 = select.z850
  iz = select.iz
  izerr = select.izerr

  int23s = findgen(121)*0.005+0.6
  scatter = int23s*0.
  ;;find a 0.2 mag wide window minimizing the measured scatter
  for i=0, n_elements(int23s)-1 do begin
     model = int23s[i] + slope*(z850-23)
     resid = iz - model
     wgals = where(abs(resid) lt 0.1)
     if wgals[0] eq -1 then scatter[i]=10. else $
        if n_elements(wgals) lt 7 then scatter[i]=10. else $
           scatter[i] = biweight_scale( resid[wgals], /zero )
  endfor

  smin = min( scatter )
  index = where( scatter eq smin ) ;;could be more than one!
  int23 = mean(int23s[index])

  ;;now use galaxies in this window to fit a red sequence
  model = int23 + slope*(z850-23)
  resid = iz - model
  wgood = where(abs(resid) lt 0.15)
  resid = resid[wgood]
  err = izerr[wgood]

;fit using a linear least squares
;  correction = total(resid/(err^2+0.03^2))/total(1./(err^2+0.03^2))
;fit using the biweight
  correction = biweight_location(resid)
  int23 += correction

  ;;now use this fit to estimate a sigma
  model = int23 + slope*(z850-23)
  resid = iz - model
  wgood = where(abs(resid) lt 0.15)
  resid = resid[wgood]
  err = izerr[wgood]
  mscatter = biweight_scale( resid, /zero, weight=weight )
  inliers = where( weight ge 0.1, complement=outliers )
;  mscatter = biweight_scale( resid[inliers], /zero )
  mmscatter = errscat2(izerr[wgood[inliers]])
  iscatter = sqrt(mscatter^2-mmscatter^2)

  whist = where( gals.select $
                 and gals.z850 ge magrange[0] $
                 and gals.z850 le magrange[1] )
  izhist = gals[whist].iz
  z850hist = gals[whist].z850
  cresid = izhist - (int23+slope*(z850hist-23))
  plothist, cresid, xhist, yhist, /noplot, xrange=[-2.,2.], bin=0.05, /nan
  histmax = max(yhist)
  plothist, cresid, xhist, yhist, peak=1, /noplot, xrange=[-2.,2.], bin=0.05, /nan

  gals.rsresid = gals.iz - int23 - slope*(gals.z850-23.)
  @define_structs2
  fit = RSfit0
  fit.slope=slope
  fit.intercept = int23
  fit.mscatter = mscatter
  fit.scatter_corr = mmscatter
  fit.iscatter=iscatter
  fit.ngals=n_elements(resid)
  fit.xhist=ptr_new(xhist)
  fit.yhist=ptr_new(yhist)
  fit.histmax=histmax
  if inliers[0] ne -1 then fit.inliers=ptr_new(wselect[wgood[inliers]]) $
  else fit.inliers=ptr_new()
  if outliers[0] ne -1 then fit.outliers=ptr_new(wselect[wgood[outliers]]) $
  else fit.outliers=ptr_new()
  fit.magrange=magrange
  fit.bluelimit=bluelimit
  cmd.fit = ptr_new(fit)
  ptr_free, cmd.gals
  cmd.gals = ptr_new(gals)

  if n_elements(diagnostic) ne 0 then begin
     thisDevice = !d.name
     set_plot, 'ps'
     device, /encapsulated, xsize=5, ysize=3, /inches
     device, filename=diagnostic+'scatter.eps'
     plot, int23s, scatter, thick=2, xtitle='int23s', ytitle='scatter', $
           yrange=[0., 0.2], ystyle=1, xthick=2, ythick=2, xrange=[0.6,1.2]
     oplot, [int23], [mscatter], ps=symcat(16), symsize=1.5
     device, /close
     set_plot, thisDevice
  endif
  


end
