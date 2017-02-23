Function reportscatter, $
   int23s, $
   slope, $
   iz, $
   z850, $
   weights=weights, $
   niter=niter, $
   c=c

  if n_elements(niter) eq 0 then niter=4
  if n_elements(c) eq 0 then c=9.

  scatter=findgen(n_elements(int23s))
  weights=findgen(n_elements(int23s), n_elements(z850))
  for i=0, n_elements(int23s)-1 do begin
     model = int23s[i] + slope * ( z850 - 23 )
     resid = iz - model
     scatter[i] = biweight_scale( resid, /zero, MAD=0.06/0.6745, $
                                  weight=weight, niter=niter, c=c )
     weights[i,*]=weight
  endfor
  return, scatter
end

Pro fit_rs8, $
   cmd, $
   slope, $
   int23guess, $
   magrange=magrange, $
   diagnostic=diagnostic

  if n_elements(magrange) eq 0 then magrange=[19.,25.]

  gals=*cmd.gals
  resid0 = (int23guess+slope*(gals.z850 - 23.))-gals.iz
  ;; pick morph E-types in appropriate range
  wselect = where( gals.select $
                   and gals.z850 ge magrange[0] $
                   and gals.z850 le magrange[1] $
                   and abs(resid0) le 0.25 )
  select = gals[wselect]
  z850 = select.z850
  iz = select.iz
  izerr = select.izerr

  ;; test intercepts
  int23s = (findgen(101)-50)*0.005+int23guess
  scatter = reportscatter( int23s, slope, iz, z850, $
                            weights=weights, niter=1, c=6 )
  smin = min( scatter, /nan ) ;; find best scatter
  index = where( scatter eq smin )  ;; more than one minimum?
  index = index[n_elements(index)-1]
  ;; in that rare case, choose the reddest possibility
  int23 = int23s[index]
  resid = iz-(int23+slope*(z850-23))
  wgood = where(abs(resid) lt 0.15)
  mscatter = reportscatter( int23, slope, iz[wgood], z850[wgood], $
                            niter=1, weights=weights )
  inliers = where( weights ne 0., complement=outliers )
  scatter_corr=errscat2(izerr[wgood[inliers]])

  whist = where( gals.select $
                 and gals.z850 ge magrange[0] $
                 and gals.z850 le magrange[1] )
  izhist = gals[whist].iz
  z850hist = gals[whist].z850
  cresid = izhist - (int23+slope*(z850hist-23))
  plothist, cresid, xhist, yhist, /noplot, xrange=[-2.,2.], bin=0.05, /nan
  histmax = max(yhist)
  plothist, cresid, xhist, yhist, peak=1, /noplot, xrange=[-2.,2.], bin=0.05, /nan
  gals.rsresid=gals.iz-int23-slope*(gals.z850-23.)

  ;;estimate errors in scatter
  ;; make KDE for residuals first
  xs = (findgen(401)/400)*0.4-0.2
  ys = xs*0.
  for i=0, n_elements(wgood)-1 do begin
     err=gals[wselect[wgood[i]]].izerr
     resid=gals[wselect[wgood[i]]].rsresid
     ys += 1./sqrt(2*!dpi*err^2)*exp(-(xs-resid)^2./(2*err^2))
  endfor
  F=total(ys,/cumulative)
  F /= max(F)
  mscat = findgen(1000)
  for i=0, 999 do begin
     counter, i+1, 1000
     u=randomu( seed, n_elements(wgood) )
     r = interpol( xs, F, u ) ;;r are distributed like KDE calculated above
     mscat[i] = biweight_scale( r, /zero, MAD=0.06/0.6745, $
                                niter=1, weight=weight )
  endfor
  iscat2 = mscat^2-scatter_corr^2
  w=where(iscat2 ge 0)
  mscat = mscat[w]
  

  @define_structs2
  fit = RSfit0
  fit.slope=slope
  fit.intercept = int23
  fit.mscatter = mscatter
  fit.mscatter_err = stdev(mscat)
  fit.scatter_corr = scatter_corr
  fit.iscatter = sqrt(mscatter^2 - scatter_corr^2)
  fit.iscatter_err = stdev(sqrt(mscat^2-scatter_corr^2))
  fit.ngals = n_elements(inliers)
  fit.xhist=ptr_new(xhist)
  fit.yhist=ptr_new(yhist)
  fit.histmax=histmax
  if inliers[0] ne -1 then fit.inliers=ptr_new(wselect[wgood[inliers]]) $
  else fit.inliers=ptr_new()
  if outliers[0] ne -1 then fit.outliers=ptr_new(wselect[wgood[outliers]]) $
  else fit.outliers=ptr_new()
  fit.magrange=magrange
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
     colora = iz - slope*(z850-23.)
     int23a = biweight_location(colora)
     residb = iz-(int23a+slope*(z850-23.))
     wgoodb = where(abs(residb) lt 0.15)
     mscata = reportscatter( int23a, slope, iz[wgoodb], z850[wgoodb], $
                             niter=1, weights=weights )
     oplot, [int23a], [mscata], ps=symcat(16), symsize=1.5, color=fsc_color('blue',254)

     device, /close
     set_plot, thisDevice
  endif
  
end
