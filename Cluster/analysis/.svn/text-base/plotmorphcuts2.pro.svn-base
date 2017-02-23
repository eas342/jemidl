Pro plotmorphcuts2, objs, hardcopy=hardcopy
  if n_elements(hardcopy) ne 0 then begin
     thisDevice = !d.name
     set_plot, 'ps'
     device, /encapsulated, /color, bits_per_pixel=8, landscape=0
     device, xsize=4.5, ysize=4.5, /inches
     device, filename=hardcopy
     thick=4
     charthick=4
     xthick=4
     ythick=4
  endif else begin
     window, 0, xsize=450, ysize=450
     thick=1
     charthick=1
     xthick=1
     ythick=1
  endelse
  resolve_obj, objs[0]
  colors = Obj_New("IDLgrPalette")
  colors->LoadCT, 3
  colors->GetProperty, Red=r, Green=g, Blue=b
  Obj_Destroy, colors
  TVLCT, r, g, b
  red = fsc_color('red', !d.table_size-2)
  grey = fsc_color('light gray', !d.table_size-3)
  blue = fsc_color('blue', !d.table_size-4)
  green = fsc_color('green', !d.table_size-5)
  
  clusters = ['A'  , 'B' , 'C' , 'D' , 'E' , 'F' , 'H' , 'K' , 'N' , 'R' , 'T' , 'U'  , 'V' , 'X' , 'Y' , 'Z'  ]
  zmin =     [1.44 , 1.10, 0.95, 1.00, 1.01, 1.10, 1.22, 1.39, 1.01, 1.2 , 0.96, 1.02 , 0.89, 1.09, 1.22, 1.375]
  zmax =     [1.476, 1.14, 0.99, 1.03, 1.04, 1.13, 1.25, 1.42, 1.03, 1.23, 0.98, 1.055, 0.92, 1.12, 1.25, 1.40 ]

  clusterids=strarr(n_elements(objs))
  for icluster=0, n_elements(objs)-1 do begin
     clusterids[icluster]=objs[icluster]->extract('clusterid')
  endfor

  for icluster=0, n_elements(clusters)-1 do begin
     w=(where(clusters[icluster] eq clusterids))[0]
     if w[0] eq -1 then continue
     s1=objs[w]->summary()
     w1=where(s1.z gt zmin[icluster] and s1.z lt zmax[icluster])
     if w1[0] eq -1 then continue
     if n_elements(s) eq 0 then s=s1[w1] $
     else s=[s,s1[w1]]
  endfor
  wno2 = where(strpos(s.comment,'oII') eq -1 and s.comment ne '?' and s.comment ne '??')
  wo2  = where(strpos(s.comment,'oII') ne -1 and s.comment ne '?' and s.comment ne '??')
  readcol, '/home/scpdata02/clusters/supernova.txt', $
           SNname, SNnickname, SNz, SNhostra, SNhostdec, $
           format='A,A,X,X,X,X,D,A,A', /silent
  wnohost=where(SNhostra eq '??', complement=whost)
  SNname=SNname[whost]
  SNnickname=SNnickname[whost]
  SNz=SNz[whost]
  SNhostra=SNhostra[whost]
  SNhostdec=SNhostdec[whost]

  ismem = member(SNname, $
                 'CL-'+['H-005','K-000','K-018','D-000','R-012','U-004' $
                        ,'B-004','V-006','F-012'])
  wsne = where(ismem)

  plotsym, 0, color=red, thick=thick

  plot, s[wno2].gini2, s[wno2].asym2, ps=8, xrange=[0.1, 0.7], yrange=[-0.05, 0.5], $
        xtitle='Gini Coefficient', ytitle='Asymmetry', $
        xstyle=1, ystyle=1, xthick=xthick, ythick=ythick, charthick=charthick, $
        title='Comparison of [OII] to automated morphology.'
  plotsym, 0, color=blue, thick=thick
  oplot, s[wo2].gini2, s[wo2].asym2, ps=8
  for iwsne=0, n_elements(wsne)-1 do begin
     get_coords, coords, instring=SNhostra[wsne[iwsne]]+' '+SNhostdec[wsne[iwsne]]
     alpha=coords[0]*15.
     delta=coords[1]
     dist2 = (s.alphawin_j2000 - alpha)^2*cos(delta*!dpi/180.d)^2 $
             + (s.deltawin_j2000 - delta)^2
     md2 = min( dist2, m)*3600.d^2
     if member(m, wno2) then begin
        plotsym, 0, color=red, /fill, thick=thick
     endif else begin
        plotsym, 0, color=blue, /fill, thick=thick
     endelse
     oplot, [s[m].gini2], [s[m].asym2], ps=8
     xyouts, [s[m].gini2], [s[m].asym2], strmid(SNname[wsne[iwsne]],3), charthick=charthick
     print, strmid(SNname[wsne[iwsne]],3)+' = '+s[m].clusterid+'.'+strtrim(s[m].galid)+' '+s[m].comment
     print, s[m].asym2, s[m].gini2
  endfor
  oplot, [0.44, 0.44], [-0.05, 0.08]
  oplot, [0.44, 0.7],  [0.08, 0.08]


  if n_elements(hardcopy) ne 0 then begin
     device, /close_file
     set_plot, thisDevice
  endif

end
