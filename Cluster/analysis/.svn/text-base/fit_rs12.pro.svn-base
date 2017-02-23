Pro fit_rs12, $
   clustercmd, $
   bkgcmd, $
   clusterarea, $
   bkgarea, $
   slope, $
   int23, $
   rsresidrange=rsresidrange, $
   magrange=magrange, $
   diagnostic=diagnostic, $
   start=start, $
   scale=scale

  ;;select galaxies to fit and get residuals
  clustergals = *clustercmd.gals
  clusterrsresid = clustergals.iz - slope*(clustergals.z850-23.)
  wcluster = where( clustergals.select $
                    and clustergals.z850 ge magrange[0] $
                    and clustergals.z850 le magrange[1] $
                    and clusterrsresid ge rsresidrange[0] $
                    and clusterrsresid le rsresidrange[1] )
  clusterdata = { rsresid:clusterrsresid[wcluster], $
                  izerr:clustergals[wcluster].izerr }
  correction = errscat2(clusterdata.izerr)

  ;;select background galaxies to fit and get residuals
  bkggals = *bkgcmd.gals
  bkgrsresid = bkggals.iz - slope*(bkggals.z850-23.)
  wbkg = where( bkggals.select $
                and bkggals.z850 ge magrange[0] $
                and bkggals.z850 le magrange[1] $
                and bkgrsresid ge rsresidrange[0] $
                and bkgrsresid le rsresidrange[1] )
  bkgdata = { rsresid:bkgrsresid[wbkg], $
              izerr:bkggals[wbkg].izerr }
  ;;initial values for fit
  scatter = 0.05
  ngal = n_elements(clusterdata.rsresid)-n_elements(bkgdata.rsresid)/bkgarea*clusterarea > 1.
  bkgconst = 2.0
  bkgslope = -0.33
  if n_elements(start) eq 0 then start = [int23, scatter, ngal, bkgconst, bkgslope, -15., -50.]
  if n_elements(scale) eq 0 then scale = [0.01, 0.002, ngal*0.1 > 1., 0.1, 0.1, 5., 10.]
  
  rslike = obj_new( 'rslike', $
                    clusterdata = ptr_new(clusterdata), $
                    bkgdata = ptr_new(bkgdata), $
                    izerr = errscat2(clusterdata.izerr), $
                    clusterarea = clusterarea, $
                    bkgarea = bkgarea, $
                    rsresidrange = rsresidrange )
  result = start
  for i=0, 9 do begin
     result1 = metropolis( 'LogLikelihood', result, scale, object=rslike, miniter=1000L, fval=fval )
;     print, fval
     result = downhillsimplex( 'LogLikelihood', result1, scale/3., object=rslike, $
                               maxiter=1000L, miniter = 100L, fval=fval )
;     print, fval[0]
  endfor
  print, result

  ;; error in scatter measurement
  scattry = findgen(1000)*0.0001
  like = fltarr(1000)
  for i=0, 999 do like[i] = rslike->LogLikelihood( [result[0], scattry[i], result[2:*]] )
  whigh = where(scattry ge result[1])
  wlow = where(scattry le result[1])
  iscatUL = 0.2
  iscatLL = -0.1
  if n_elements(whigh) gt 1 and result[2] gt 0. $
  then iscatUL = (interpol(scattry[whigh], like[whigh], min(fval)+0.5)) < 0.2
  if n_elements(wlow) gt 1 and result[2] gt 0. $
  then iscatLL = (interpol(scattry[wlow], like[wlow], min(fval)+0.5)) > (-0.1)
  
  if ~finite(iscatLL) then iscatLL = -0.1
  if ~finite(iscatUL) then iscatLL = 0.2

  junk = rslike->LogLikelihood(result)

  ;;histograms
  whist = where( clustergals.select $
                 and finite(clusterrsresid) $
                 and clustergals.z850 ge magrange[0] $
                 and clustergals.z850 le magrange[1] )
  izhist = clustergals[whist].iz
  z850hist = clustergals[whist].z850
  cresid = izhist - (result[0]+slope*(z850hist-23.))
  plothist, cresid, xhist, yhist, /noplot, xrange=[-2.,2.], bin=0.05, /halfbin
  histmax = max(yhist)
  plothist, cresid, xhist, yhist, peak=1, /noplot, xrange=[-2., 2.], bin=0.05, /halfbin
  clustergals.rsresid=clustergals.iz-result[0]-slope*(clustergals.z850-23.)
  ;;background histogram
  wplotbkg = where( bkggals.select $
                    and finite(bkgrsresid) $
                    and bkggals.z850 ge magrange[0] $
                    and bkggals.z850 le magrange[1] )
  bresid = (bkggals.iz - (result[0]+slope*(bkggals.z850-23.)))[wplotbkg]
  plothist, bresid, bin=0.05, bkgxhist, bkgyhist, /noplot, /halfbin
  bkgyhist *= clusterarea/bkgarea/histmax
  ;;cluster fit
  xplot = findgen(500)/499.*(rsresidrange[1]-rsresidrange[0])+rsresidrange[0]
  yplot = rslike->ClusterIntegrand(xplot)*0.05
  ;;bkg fit
  bkgxplot = xplot
  bkgyplot = rslike->BkgIntegrand(xplot)*0.05*clusterarea/bkgarea
  
  winliers = where( abs((clusterdata.rsresid-result[0])) $
                    lt 2.*sqrt(result[1]^2+clusterdata.izerr^2),  $
                    complement = woutliers )
  if winliers[0] ne -1 then inliers = wcluster[winliers]
  if woutliers[0] ne -1 then outliers = wcluster[woutliers]



  @define_structs2
  fit = RSfit0
  fit.slope=slope
  fit.intercept = result[0]
  fit.mscatter = sqrt(result[1]^2+correction^2)
  fit.mscatter_err = 0.
  fit.scatter_corr = correction
  fit.iscatter = abs(result[1])
  fit.iscatter_err = 0.
  fit.iscatLL = iscatLL
  fit.iscatUL = iscatUL
  fit.ngals = result[2]
  fit.xhist=ptr_new(xhist)
  fit.yhist=ptr_new(yhist)
  fit.histmax=histmax
  fit.bkgxhist = ptr_new(bkgxhist)
  fit.bkgyhist = ptr_new(bkgyhist)
  fit.xplot = ptr_new(xplot)
  fit.yplot = ptr_new(yplot)
  fit.bkgxplot = ptr_new(bkgxplot)
  fit.bkgyplot = ptr_new(bkgyplot)
  if n_elements(inliers) ne 0 then fit.inliers=ptr_new(inliers) $
  else fit.inliers=ptr_new()
  if n_elements(outliers) ne 0 then fit.outliers=ptr_new(outliers) $
  else fit.outliers=ptr_new()
  fit.magrange=magrange
  clustercmd.fit = ptr_new(fit)
  ptr_free, clustercmd.gals
  clustercmd.gals = ptr_new(clustergals)

  if n_elements(diagnostic) ne 0 then begin
     thisDevice = !d.name
     set_plot, 'ps'
     device, /encapsulated, xsize=5, ysize=3, /inches
     device, filename=diagnostic
     
     plot, scattry, 2*(like-min(like)), xtitle='Intrinsic scatter', ytitle='2 * Log Likelihood', $
           thick=3, charthick=3, xthick=3, ythick=3, yrange=[0,5], ystyle=1, title=diagnostic
     oplot, scattry, scattry*0+1, thick=3

     device, /close
     set_plot, thisDevice
  endif
  
end
