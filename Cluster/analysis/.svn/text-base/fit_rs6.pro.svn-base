Pro fit_rs6, cmd, slope, magrange=magrange, bluelimit=bluelimit
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

  int23s = findgen(141)*0.005+0.6
  ngals = int23s*0.
  ;;find a 0.2 mag wide window maximizing ngal
  for i=0, n_elements(int23s)-1 do begin
     model = int23s[i] + slope*(z850-23)
     resid = iz - model
     ngals[i] = total(abs(resid) lt 0.1)
  endfor
  nmax = max( ngals )
  index = where( ngals eq nmax ) ;;could be more than one!
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
  inliers = where( weight ge 0.8, complement=outliers )
  mmscatter = errscat2(izerr[wgood[inliers]])
  iscatter = sqrt(mscatter^2-mmscatter^2)

  izall = gals[where(gals.select)].iz
  z850all = gals[where(gals.select)].z850
  c = izall - slope*(z850all-25)
  plothist, c, xhist, yhist, peak=1, /noplot, xrange=[-0.5,1.5], bin=0.05, /nan

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
  if inliers[0] ne -1 then fit.inliers=ptr_new(wselect[wgood[inliers]]) $
  else fit.inliers=ptr_new()
  if outliers[0] ne -1 then fit.outliers=ptr_new(wselect[wgood[outliers]]) $
  else fit.outliers=ptr_new()
  fit.magrange=magrange
  fit.bluelimit=bluelimit
  cmd.fit = ptr_new(fit)
  ptr_free, cmd.gals
  cmd.gals = ptr_new(gals)
end
