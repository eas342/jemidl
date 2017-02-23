Pro rsmetropoliswrapper, obj, start=start, scale=scale, magrange=magrange, colorrange=colorrange
  common JEM$_rssimplexwrapper, bkgdata
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; int23, slope,  sig,  phi*,  m*, alpha,  d,  e,  f,    g
;  if n_elements(start) eq 0 then start = [0.946, -0.028, 0.028, 5., 22., -1., -7., 0.7, -3., -0.007]
  if n_elements(start) eq 0 then start = [0.85, -0.05,  0.1, 5., 22.,   -1., -10.,  0.7, -3.7, 0.0009]
  if n_elements(scale) eq 0 then scale = [0.05,  0.01, 0.02, 1., 0.3, 0.03, 1.0, 0.1, 0.2, 0.002]
  if n_elements(magrange) eq 0 then magrange = [20., 24.]  ;;these defaults should be more dynamic
  if n_elements(colorrange) eq 0 then colorrange = [0.6, 1.2]

  resolve_obj, obj
  zcluster = obj->extract('zcluster')
  s=obj->summary()
  cmd=cmd( obj, /xpsf, /nofit, err=7 )
  gals=*cmd.gals
  w=where( gals.zmag gt magrange[0] $
           and gals.zmag le magrange[1] $
           and gals.iz gt colorrange[0] $
           and gals.iz le colorrange[1] $
           and gals.rad and ~gals.star )
  gals=gals[w]
  clusterdata = { mag:0., color:0., msigma:0., csigma:0. }
  clusterdata = replicate( clusterdata, n_elements(gals) )
  clusterdata.mag = gals.zmag
  clusterdata.color = gals.iz
  clusterdata.msigma = 0.
  clusterdata.csigma = gals.iz_err
  ptr_free, cmd.gals

  iexptime = mean( s[w].iexptime )
  zexptime = mean( s[w].zexptime )

  iskysigma = mean( s[w].isky_sigma )
  zskysigma = mean( s[w].zsky_sigma )

  
  if n_elements(bkgdata) eq 0 then begin
     restore, '/home/scpdata03/goods/cmd/cmds.sav'
     bkggals = *(*cmds[0]).gals
     for i=1, n_elements(cmds)-1 do begin
        bkggals=[bkggals, *(*cmds[i]).gals]
     endfor
     w=where( bkggals.zmag gt magrange[0] $
              and bkggals.zmag le magrange[1] $
              and bkggals.iz gt colorrange[0] $
              and bkggals.iz le colorrange[1] $
              and ~bkggals.star )
     bkggals = bkggals[w]
     bkgdata = {mag:0., color:0.}
     bkgdata = replicate( bkgdata, n_elements(bkggals) )
     bkgdata.mag = bkggals.zmag
     bkgdata.color = bkggals.iz
     for i=0, n_elements(cmds)-1 do begin
        ptr_free, (*cmds[i]).gals
     endfor
     ptr_free, cmds
  endif

  rMpc = 0.8
  rarcmin = rMpc / (lumdist( zcluster, /silent )/(1.+zcluster)^2) * 180./!pi * 60.
  clusterarea = rarcmin^2*!pi
  bkgarea = 1.4^2*!pi*30

  print
  print, format='(%"Cluster galaxies:            %4i")', n_elements(clusterdata)
  print, format='(%"Background galaxies:         %4i")', n_elements(bkgdata)
  print, format='(%"Cluster galaxy density:     %5.2f")', n_elements(clusterdata)/clusterarea
  print, format='(%"Background galaxy density:  %5.2f")', n_elements(bkgdata)/bkgarea
  print, format='(%"Cluster area:             %7.2f")', clusterarea
  print, format='(%"Background area:          %7.2f")', bkgarea
  print, format='(%"Cluster excess:           %7.2f")', n_elements(clusterdata)-n_elements(bkgdata)/bkgarea*clusterarea
  print

  rslike = obj_new( 'rslike', $
                    clusterdata=ptr_new(clusterdata), $
                    bkgdata=ptr_new(bkgdata), $
                    clusterarea=clusterarea, $
                    bkgarea=bkgarea, $
                    magrange=magrange, $
                    colorrange=colorrange, $
                    iexptime=iexptime, $
                    zexptime=zexptime, $
                    iskysigma=iskysigma, $
                    zskysigma=zskysigma )

  junk = rslike->loglikelihood(start)

  result = metropolis( 'LogLikelihood', start, scale, object = rslike, miniter=100, /verbose )
  for i=0, 10 do begin
     result = metropolis( 'LogLikelihood', result, scale, object = rslike, $
                               miniter=100, fval=fval, /verbose )

     print
     print
     print, "Results!"
     print, format='(%"CMRint23:  %8.4f (%8.4f)")', result[0], start[0]
     print, format='(%"CMRslope:  %8.4f (%8.4f)")', result[1], start[1]
     print, format='(%"CMRsigma:  %8.4f (%8.4f)")', result[2], start[2]
     print, format='(%"phi_star:  %8.4f (%8.4f)")', result[3], start[3]
     print, format='(%"m_star:    %8.4f (%8.4f)")', result[4], start[4]
     print, format='(%"alpha:     %8.4f (%8.4f)")', result[5], start[5]
     print, format='(%"bkg_d:     %8.4f (%8.4f)")', result[6], start[6]
     print, format='(%"bkg_e:     %8.4f (%8.4f)")', result[7], start[7]
     print, format='(%"bkg_f:     %8.4f (%8.4f)")', result[8], start[8]
     print, format='(%"bkg_g:     %8.4f (%8.4f)")', result[9], start[9]
     print, format='(%"loglikelihood:  %10.2f")', fval[0]
     print
  endfor

  obj_destroy, rslike
end
