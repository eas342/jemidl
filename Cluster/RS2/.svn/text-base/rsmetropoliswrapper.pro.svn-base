Function rsmetropoliswrapper, $
   obj, $
   slope, $
   start=start, $
   scale=scale, $
   magrange=magrange, $
   rsresidrange=rsresidrange, $
   plothist=plothist
  
  common JEM$_rssimplexwrapper, bkgdata, bkgarea, bkggals, bkgrsresid
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; int23,  sig,   N,   d,    e
  if n_elements(start) eq 0 then start = [0.95,   0.07,  10., 0.2,  0.]
  if n_elements(scale) eq 0 then scale = [0.02,   0.01,  1.,  0.01, 0.001]
  if n_elements(rsresidrange) eq 0 then rsresidrange=[0.6,1.2]
  if n_elements(magrange) eq 0 then magrange=[20.,24]

  resolve_obj, obj
  zcluster = obj->extract('zcluster')
  clustercmd=cmd2( obj, /xpsf )
  selectgals, clustercmd
  clustergals=*clustercmd.gals
  clusterrsresid = clustergals.iz - slope*(clustergals.z850-23.)

  wcluster =  where( clustergals.select $
                     and clustergals.z850 ge magrange[0] $
                     and clustergals.z850 le magrange[1] $
                     and clusterrsresid ge rsresidrange[0] $
                     and clusterrsresid le rsresidrange[1] )
  clusterdata = { rsresid:clusterrsresid[wcluster], $
                  izerr:clustergals[wcluster].izerr }
  rMpc = 0.65
  rarcmin = rMpc / (lumdist( zcluster, /silent )/(1.+zcluster)^2) * 180./!pi * 60.
  clusterarea = rarcmin^2*!pi

  ;;Load GOODS data if not explicitly passed
  if n_elements(bkgdata) eq 0 then begin
     restore, '/home/scpdata03/goods/cmd/cmds.sav'
     bkggals = *(*cmds[0]).gals
     for i=1, 29 do begin
        cmd1 = *cmds[i]
        bkggals = [bkggals, *cmd1.gals]
     endfor
     bkgcmd = {clusterid:'G', clustername:'GOODS', zcluster:0.d, $
               fit:ptr_new(), gals:ptr_new(bkggals)}
     selectgals, bkgcmd
     bkgarea = 1.4^2*!pi*30.
     
     bkggals = *bkgcmd.gals
     bkgrsresid = bkggals.iz - slope*(bkggals.z850-23.)
     wbkg = where( bkggals.select $
                   and bkggals.z850 ge magrange[0] $
                   and bkggals.z850 le magrange[1] $
                   and bkgrsresid ge rsresidrange[0] $
                   and bkgrsresid le rsresidrange[1] )
     bkgdata = { rsresid:bkgrsresid[wbkg], $
                 izerr:bkggals[wbkg].izerr }
     
     for i=0, n_elements(cmds)-1 do begin
        ptr_free, (*cmds[i]).gals
     endfor
     ptr_free, cmds
  endif
  
  print
  print, format='(%"Cluster galaxies:            %4i")', $
         n_elements(clusterdata.rsresid)
  print, format='(%"Background galaxies:         %4i")', $
         n_elements(bkgdata.rsresid)
  print, format='(%"Cluster galaxy density:    %6.3f")', $
         n_elements(clusterdata.rsresid)/clusterarea
  print, format='(%"Background galaxy density: %6.3f")', $
         n_elements(bkgdata.rsresid)/bkgarea
  print, format='(%"Cluster excess:           %7.2f")', $
         n_elements(clusterdata.rsresid)-n_elements(bkgdata.rsresid)/bkgarea*clusterarea
  print

  rslike = obj_new( 'rslike', $
                    clusterdata=ptr_new(clusterdata), $
                    bkgdata=ptr_new(bkgdata), $
                    clusterarea=clusterarea, $
                    bkgarea=bkgarea, $
                    rsresidrange=rsresidrange )

  junk = rslike->loglikelihood(start)

  start[2] = n_elements(clusterdata.rsresid)-n_elements(bkgdata.rsresid)/bkgarea*clusterarea > 1.
  start[2] /= clusterarea

  result = metropolis( 'LogLikelihood', start, scale, object = rslike, miniter=10000, fval=fval )

  print
  print
  print, "Results!"
  print, format='(%"CMRint23:  %8.4f (%8.4f)")', result[0], start[0]
  print, format='(%"CMRsigma:  %8.4f (%8.4f)")', result[1], start[1]
  print, format='(%"N:         %8.4f (%8.4f)")', result[2]*clusterarea, start[2]*clusterarea
  print, format='(%"bkg_d:     %8.4f (%8.4f)")', result[3], start[3]
  print, format='(%"bkg_e:     %8.4f (%8.4f)")', result[4], start[4]
  print, format='(%"loglikelihood:  %10.2f")', fval[0]
  print
  junk = rslike->loglikelihood(result, /verbose)

  s=findgen(1000)*0.0001
  ss=fltarr(1000)
  for i=0, 999 do begin
     ss[i] = rslike->loglikelihood([result[0], s[i], result[2:*]])
  endfor
  window, 1, xsize=800, ysize=500
  plot, s, ss-min(ss)
  oplot, s, s*0+1
  junk = rslike->loglikelihood(result)

  if keyword_set(plothist) then begin
     window, 0, xsize=800, ysize=500
     junk = fsc_color(/all, color=ctable)
     wplotcluster = where( clustergals.z850 ge magrange[0] $
                           and clustergals.z850 le magrange[1] $
                           and finite(clusterrsresid) $
                           and clustergals.select )
     cresid = clustergals.iz - result[0] - slope*(clustergals.z850-23.)
     cresid = cresid[wplotcluster]
     plothist, cresid, x, y, bin=0.05, /halfbin, /noplot
     center = result[0]+slope*(25-23.)
     plot, [0], /nodata, xrange=[-0.5,1.5], yrange=[0,max(y)+5]
     oplot, x+center, y, ps=10, color=ctable.red
;     plothist, cresid, bin=0.05, color=ctable.yellow, /halfbin, $
;               xrange=[-0.5, 1.5], yrange=[0,max(y)+5], /over
;     plotslopehist, x-0.025, y, 0., color=ctable.green, xrange=[-0.5,1.5], yrange=[0,max(y)+5]
     wplotbkg = where( bkggals.z850 ge magrange[0] $
                       and bkggals.z850 le magrange[1] $
                       and finite(bkgrsresid) $
                       and bkggals.select )
     plothist, bkgrsresid[wplotbkg], bin=0.05, x, y, /noplot, /halfbin
     oplot, x, y*clusterarea/bkgarea, ps=10, color=ctable.blue
     xs = findgen(500)/499.*(rsresidrange[1]-rsresidrange[0])+rsresidrange[0]
     yplot = rslike->ClusterIntegrand(xs)*0.05
     oplot, xs, yplot, color=ctable.magenta
     yplot = rslike->BkgIntegrand(xs)*0.05*clusterarea/bkgarea
     oplot, xs, yplot, color=ctable.cyan
  endif

  obj_destroy, rslike

  return, result
end
