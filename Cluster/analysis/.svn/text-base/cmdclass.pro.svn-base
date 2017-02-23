Function cmdclass, obj, colorrange=colorrange, magrange=magrange, $
                   zrange=zrange, slope=slope, morph=morph, $
                   clusterrad=clusterrad, plotrange=plotrange, $
                   _extra=extra
                   
  if n_elements(clusterrad) eq 0 then clusterrad = 0.8
  if n_elements(morph) eq 0 then morph=4
  if n_elements(slope) eq 0 then slope=0.
  if n_elements(colorrange) eq 0 then colorrange=[0.8, 1.1]
  if n_elements(magrange) eq 0 then magrange=[20.,23.5]
  
  if keyword_set(ps) then begin
     thick=3
  endif else begin
     thick=1
  endelse

  if keyword_set(plotrange) then begin
     junk = fsc_color(/allcolors, colorstructure=ctable)
     oplot, [19,25], ([19,25]-23)*slope+colorrange[0], color=ctable.blue, linestyle=2, thick=thick
     oplot, [19,25], ([19,25]-23)*slope+colorrange[1], color=ctable.blue, linestyle=2, thick=thick
     oplot, [1.,1.]*magrange[0], [-0.5, 1.5], color=ctable.blue, linestyle=2, thick=thick
     oplot, [1.,1.]*magrange[1], [-0.5, 1.5], color=ctable.blue, linestyle=2, thick=thick
  endif

  
  resolve_obj, obj
  s=obj->summary()
  c=obj->color(_extra=extra, /silent)
  m=s.zmag_auto

  ;; MW dust correction
  ebv = obj->extract('ebv')
  A_i = 1.973*ebv
  A_z = 1.472*ebv
  m -= A_z
  c -= (A_i - A_z)

  ;; SExtractor MAG_AUTO bias correction
  ;m -= (m-21)*0.1/3.+0.2
  m += m*(-0.0422781)+0.635273


  ;; clustercentric radius cut
  zcluster=obj->extract('zcluster')
  ra_center=obj->extract('ra')
  dec_center=obj->extract('dec')
  ra=s.alphawin_j2000
  dec=s.deltawin_j2000
  radius2 = (ra_center-ra)^2*(cos(dec_center*!dpi/180.d))^2 $
            + (dec_center-dec)^2 ;; clustercentric distance^2 in degrees^2
  rthresh = clusterrad/(lumdist( zcluster, /silent )/(1.+zcluster)^2) ;; in rad
  rthresh *= 180.d/!dpi ;; rad->degrees
  radcheck = radius2 le rthresh^2

  ;; star cut
  starcheck = s.zflux_radius lt 2.2

  ;; morphcut
  case morph of
     1: morphcheck = s.conc ge 0.35+4.65*s.asym
     2: morphcheck = s.asym le 0.1 and s.gini ge 0.44
     3: morphcheck = s.asym le 0.08 and s.gini ge 0.44
     4: morphcheck = s.asym le 0.1 and s.gini ge 0.43
  endcase
  
  rsregioncheck = m ge magrange[0] and m le magrange[1] $
                  and c ge (colorrange[0] + (m-23)*slope) $
                  and c le (colorrange[1] + (m-23)*slope)
  
  if n_elements(zrange) eq 0 then zrange = zcluster + [-0.02, 0.02]
  specmemcheck = s.z ge zrange[0] and s.z le zrange[1]
  specnoncheck = s.z ge 0. and ~specmemcheck
  no2check = strpos( s.comment, 'oII' ) eq -1 and (specmemcheck or specnoncheck) $
             and s.comment ne '?' and s.comment ne '??'
  o2check = strpos( s.comment, 'oII' ) ne -1
  nocluecheck = s.comment eq '?' or s.comment eq '??'
  cerr = geterrs(obj, _extra=extra)
  merr = s.zmagerr_auto
  return, { rad:radcheck, $
            star:starcheck, $
            morph:morphcheck, $
            rsregion:rsregioncheck, $
            specmem:specmemcheck, $
            specnon:specnoncheck, $
            no2:no2check, $
            o2:o2check, $
            noclue:nocluecheck, $
            m:m, $
            c:c, $
            merr:merr, $
            cerr:cerr } 
end
