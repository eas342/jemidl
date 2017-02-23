Pro plotmorphcuts, objs, hardcopy=hardcopy
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
  
  clusters = ['A' ,'B'  ,'C'  ,'D'  ,'E'  ,'F'  ,'H'   ,'K'  ,'N'  ,'R'  ,'T'  ,'U'   ,'V'  ,'X'  ,'Y'  ,'Z']
  zmin =     [1.43, 1.1 , 0.94, 1.00, 1.  , 1.1 , 1.225, 1.4 , 1.01, 1.2 , 0.97, 1.02 , 0.9 , 1.05, 1.22, 1.39]
  zmax =     [1.47, 1.14, 1.  , 1.03, 1.05, 1.15, 1.235, 1.42, 1.03, 1.22, 0.98, 1.055, 0.91, 1.15, 1.25, 1.40]

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
  wno2 = where(strpos(s.comment,'oII') eq -1 or s.comment eq '?', complement=wo2)

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

  ;;0
  window, 0, xsize=800, ysize=800
  !p.multi=[0,2,2]
  plotsym, 0, color=red
  plot, s[wno2].gini, s[wno2].conc, ps=8, xrange=[0.1, 0.7], yrange=[0.1, 0.7], $
        xtitle='gini coefficient', ytitle='concentration'
  plotsym, 0, color=blue
  oplot, s[wo2].gini, s[wo2].conc, ps=8
  for iwsne=0, n_elements(wsne)-1 do begin
     get_coords, coords, instring=SNhostra[wsne[iwsne]]+' '+SNhostdec[wsne[iwsne]]
     alpha=coords[0]*15.
     delta=coords[1]
     dist2 = (s.alphawin_j2000 - alpha)^2*cos(delta*!dpi/180.d)^2 $
             + (s.deltawin_j2000 - delta)^2
     md2 = min( dist2, m)*3600.d^2
     if member(m, wno2) then begin
        color=red
     endif else begin
        color=blue
     endelse
     plotsym, 0, color=color, /fill
     oplot, [s[m].gini], [s[m].conc], ps=8
     xyouts, [s[m].gini], [s[m].conc], strmid(SNname[wsne[iwsne]],3)
     print, strmid(SNname[wsne[iwsne]],3)+' = '+s[m].clusterid+'.'+strtrim(s[m].galid)+' '+s[m].comment
     print, s[m].asym, s[m].conc, s[m].gini
  endfor
  
  ;;1
  plotsym, 0, color=red
  plot, s[wno2].asym, s[wno2].conc, ps=8, xrange=[-0.05, 0.15], yrange=[0.1, 0.7], $
        xtitle='asymmetry', ytitle='concentration'
  plotsym, 0, color=blue
  oplot, s[wo2].asym, s[wo2].conc, ps=8
;  oplot, [0., 0.1], [0.33, 0.63]
;  oplot, [0., 0.1], [0.3, 0.75]
;  oplot, [0., 0.07], [0.3, 0.615]
  oplot, [0., 0.07], [0.35, 0.675]
  for iwsne=0, n_elements(wsne)-1 do begin
     get_coords, coords, instring=SNhostra[wsne[iwsne]]+' '+SNhostdec[wsne[iwsne]]
     alpha=coords[0]*15.
     delta=coords[1]
     dist2 = (s.alphawin_j2000 - alpha)^2*cos(delta*!dpi/180.d)^2 $
             + (s.deltawin_j2000 - delta)^2
     md2 = min( dist2, m)*3600.d^2
     if member(m, wno2) then begin
        plotsym, 0, color=red, /fill
     endif else begin
        plotsym, 0, color=blue, /fill
     endelse
     oplot, [s[m].asym], [s[m].conc], ps=8
     xyouts, [s[m].asym], [s[m].conc], strmid(SNname[wsne[iwsne]],3)
  endfor


  ;;2
  plotsym, 0, color=red
  plot, s[wno2].gini, s[wno2].asym, ps=8, xrange=[0.1, 0.7], yrange=[-0.05, 0.15], $
        xtitle='gini coefficient', ytitle='asymmetry'
  plotsym, 0, color=blue
  oplot, s[wo2].gini, s[wo2].asym, ps=8
  !p.multi=0
  for iwsne=0, n_elements(wsne)-1 do begin
     get_coords, coords, instring=SNhostra[wsne[iwsne]]+' '+SNhostdec[wsne[iwsne]]
     alpha=coords[0]*15.
     delta=coords[1]
     dist2 = (s.alphawin_j2000 - alpha)^2*cos(delta*!dpi/180.d)^2 $
             + (s.deltawin_j2000 - delta)^2
     md2 = min( dist2, m)*3600.d^2
     if member(m, wno2) then begin
        plotsym, 0, color=red, /fill
     endif else begin
        plotsym, 0, color=blue, /fill
     endelse
     oplot, [s[m].gini], [s[m].asym], ps=8
     xyouts, [s[m].gini], [s[m].asym], strmid(SNname[wsne[iwsne]],3)
  endfor

  ;;;;;;;;;;;;;;
  ; hardcopy
  if keyword_set(hardcopy) then begin
     thisDevice=!d.name
     set_plot, 'ps'
     device, /encapsulated, color=1, bits_per_pixel=8, landscape=0
     device, xsize=8., ysize=8., /inches, xoffset=0, yoffset=0
     device, filename='ACG.eps'
     
     !p.multi = [0,2,2]
     
     ;;0
     plotsym, 0, color=red
     plot, s[wno2].gini, s[wno2].conc, ps=8, xrange=[0.1, 0.7], yrange=[0.1, 0.7], $
           xtitle='!5gini coefficient', ytitle='!5concentration'
     plotsym, 0, color=blue
     oplot, s[wo2].gini, s[wo2].conc, ps=8
     for iwsne=0, n_elements(wsne)-1 do begin
        get_coords, coords, instring=SNhostra[wsne[iwsne]]+' '+SNhostdec[wsne[iwsne]]
        alpha=coords[0]*15.
        delta=coords[1]
        dist2 = (s.alphawin_j2000 - alpha)^2*cos(delta*!dpi/180.d)^2 $
                + (s.deltawin_j2000 - delta)^2
        md2 = min( dist2, m)*3600.d^2
        if member(m, wno2) then begin
           plotsym, 0, color=red, /fill
        endif else begin
           plotsym, 0, color=blue, /fill
        endelse
        oplot, [s[m].gini], [s[m].conc], ps=8
     endfor
     
     
     ;;1
     plotsym, 0, color=red
     plot, s[wno2].asym, s[wno2].conc, ps=8, xrange=[-0.05, 0.15], yrange=[0.1, 0.7], $
           xtitle='!5asymmetry', ytitle='!5concentration'
     plotsym, 0, color=blue
     oplot, s[wo2].asym, s[wo2].conc, ps=8
;     oplot, [0., 0.1], [0.33, 0.63]
;     oplot, [0., 0.1], [0.3, 0.75]
     oplot, [0., 0.07], [0.3, 0.615]
     for iwsne=0, n_elements(wsne)-1 do begin
        get_coords, coords, instring=SNhostra[wsne[iwsne]]+' '+SNhostdec[wsne[iwsne]]
        alpha=coords[0]*15.
        delta=coords[1]
        dist2 = (s.alphawin_j2000 - alpha)^2*cos(delta*!dpi/180.d)^2 $
                + (s.deltawin_j2000 - delta)^2
        md2 = min( dist2, m)*3600.d^2
        if member(m, wno2) then begin
           plotsym, 0, color=red, /fill
        endif else begin
           plotsym, 0, color=blue, /fill
        endelse
        oplot, [s[m].asym], [s[m].conc], ps=8
     endfor
     
     
     ;;2
     plotsym, 0, color=red
     plot, s[wno2].gini, s[wno2].asym, ps=8, xrange=[0.1, 0.7], yrange=[-0.05, 0.15], $
           xtitle='!5gini coefficient', ytitle='!5asymmetry'
     plotsym, 0, color=blue
     oplot, s[wo2].gini, s[wo2].asym, ps=8
     for iwsne=0, n_elements(wsne)-1 do begin
        get_coords, coords, instring=SNhostra[wsne[iwsne]]+' '+SNhostdec[wsne[iwsne]]
        alpha=coords[0]*15.
        delta=coords[1]
        dist2 = (s.alphawin_j2000 - alpha)^2*cos(delta*!dpi/180.d)^2 $
                + (s.deltawin_j2000 - delta)^2
        md2 = min( dist2, m)*3600.d^2
        if member(m, wno2) then begin
           plotsym, 0, color=red, /fill
        endif else begin
           plotsym, 0, color=blue, /fill
        endelse
        oplot, [s[m].gini], [s[m].asym], ps=8
     endfor
     
     
     device, /close_file
     !p.multi=0
     
     set_plot, thisDevice
  endif
end
