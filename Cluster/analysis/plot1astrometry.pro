Pro plot1astrometry, obj, ps=ps, specmembers=specmembers, specnonmembers=specnonmembers, sne=sne, rsin=rsin
  if n_elements(circsize) eq 0 then circsize=1.
  resolve_obj, obj

  ;;;;;;;;;;;;;;;;;;;;;
  ; load colors

  tvlct, savered, savegreen, saveblue, /get
  if keyword_set(ps) then begin
     nomorphcolor = fsc_color('gray', !d.table_size-1)
     morphcolor = fsc_color('black', !d.table_size-2)
     background = fsc_color('white', !d.table_size-3)
     o2color = fsc_color('blue', !d.table_size-4)
     no2color = fsc_color('red', !d.table_size-5)
     noncolor = fsc_color('black', !d.table_size-6)
     SNcolor = fsc_color('orange', !d.table_size-7)
     nocluecolor = fsc_color('purple', !d.table_size-8)
     rscolor = fsc_color('red', !d.table_size-5)
  endif else begin
     nomorphcolor = fsc_color('dark gray', !d.table_size-1)
     morphcolor = fsc_color('white', !d.table_size-2)
     background = fsc_color('black', !d.table_size-3)
     o2color = fsc_color('blue', !d.table_size-4)
     no2color = fsc_color('red', !d.table_size-5)
     noncolor = fsc_color('white', !d.table_size-6)
     SNcolor = fsc_color('gold', !d.table_size-7)
     nocluecolor = fsc_color('purple', !d.table_size-8)
     rscolor = fsc_color('red', !d.table_size-5)
  endelse

  s=obj->summary()
  z=obj->extract('zcluster')
  colors = obj->color( _extra=extra )
  mags = s.zmag_best

  ebv = obj->extract('ebv')
  A_i = 1.973*ebv
  A_z = 1.472*ebv
  
  mags -= A_z

  colors -= (A_i - A_z)

  ra=s.alphawin_j2000
  dec=s.deltawin_j2000

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ; setup plot window

  xrange = reverse(minmax(ra))
  yrange = minmax(dec)

  plot, [0], [0], /nodata, xstyle=1, ystyle=1, $
        xrange=xrange, yrange=yrange, $
        title=string(format='(%"%s z=%5.2f")', obj->extract('clustername'), obj->extract('zcluster')), $
        color=morphcolor, xtitle='RA', ytitle='DEC'

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ; cut out likely [oII] emitters
  
  morphcheck = s.conc gt 0.3+4.5*s.asym
  wmorph = where( morphcheck, complement=wnomorph )
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ; plot all the gals...  will highlight next

  if wnomorph[0] ne -1 then begin
     plotsym, 0, circsize, color=nomorphcolor
     oplot, ra[wnomorph], dec[wnomorph], ps=8
  endif


  if wmorph[0] ne -1 then begin
     plotsym, 0, circsize, /fill, color=morphcolor
     oplot, ra[wmorph], dec[wmorph], ps=8
  endif

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; if red-sequence gals are given, color them in

  if n_elements(rsin) ne 0 then begin
     plotsym, 0, circsize, /fill, color=rscolor
     oplot, ra[rsin], dec[rsin], ps=8
  endif

  ;;;;;;;;;;;;;;;;
  ; highlight spectroscopic members

  if keyword_set(specmembers) or keyword_set(specnonmembers) then begin
     
     case obj->Extract('clusterid') of 
        'A' : zcheck = s.z ge 1.44  and s.z le 1.476
        'B' : zcheck = s.z ge 1.10  and s.z le 1.14
        'C' : zcheck = s.z ge 0.95  and s.z le 0.99
        'D' : zcheck = s.z ge 1.00  and s.z le 1.03
        'E' : zcheck = s.z ge 1.01  and s.z le 1.04
        'F' : zcheck = s.z ge 1.10  and s.z le 1.13
        'G' : zcheck = s.z ge 1.24  and s.z le 1.28
        'H' : zcheck = s.z ge 1.22  and s.z le 1.25
        'I' : zcheck = s.z ge 1.33  and s.z le 1.36
        'J' : zcheck = s.z ge 1.36  and s.z le 1.38
        'K' : zcheck = s.z ge 1.39  and s.z le 1.42
        'L' : zcheck = s.z ge 100.  ;;don't allow anything to pass
        'M' : zcheck = s.z ge 0.88  and s.z le 0.92
        'N' : zcheck = s.z ge 1.01  and s.z le 1.03
        'R' : zcheck = s.z ge 1.20  and s.z le 1.23
        'T' : zcheck = s.z ge 0.96  and s.z le 0.98
        'U' : zcheck = s.z ge 1.02  and s.z le 1.055
        'V' : zcheck = s.z ge 0.89  and s.z le 0.92
        'W' : zcheck = s.z ge 1.24  and s.z le 1.28
        'X' : zcheck = s.z ge 1.09  and s.z le 1.12
        'Y' : zcheck = s.z ge 1.22  and s.z le 1.25
        'Z' : zcheck = s.z ge 1.375  and s.z le 1.40
        else : zcheck = intarr(n_elements(s))*0
     endcase
     
     if keyword_set(specmembers) then begin
        if member(obj->extract('clusterid'), $
                  ['A','B','C','D','E','F','H','K','N','R','T','U','V','X','Y','Z']) then begin
           no2check = strpos( s.comment, 'oII' ) eq -1
        endif else begin
           no2check = intarr(n_elements(s))*0
        endelse

        nocluecheck = s.comment eq '?' or s.comment eq '??'
        wno2spec = where(zcheck and no2check and (1-nocluecheck))
        wo2spec = where(zcheck and (1-no2check))
        wnocluespec = where(zcheck and nocluecheck)
        if wo2spec[0] ne -1 then begin
           plotsym, 8, circsize*2, color=o2color
           oplot, ra[wo2spec], dec[wo2spec], ps=8
        endif
        if wno2spec[0] ne -1 then begin
           plotsym, 8, circsize*2, color=no2color
           oplot, ra[wno2spec], dec[wno2spec], ps=8
        endif
        if wnocluespec[0] ne -1 then begin
           plotsym, 8, circsize*2, color=nocluecolor
           oplot, ra[wnocluespec], dec[wnocluespec], ps=8
        endif
     endif
     if keyword_set(specnonmembers) then begin
        wnonspec = where((1-zcheck) and s.z ge 0.)
        if wnonspec[0] ne -1 then begin
           oplot, ra[wnonspec], dec[wnonspec], ps=7, symsize=circsize*1.8, color=noncolor
        endif
     endif
  endif

  if keyword_set(sne) then begin
     readcol, '/home/scpdata02/clusters/supernova.txt', $
              SNname, SNnickname, SNhostra, SNhostdec, $
              format='A,A,X,X,X,X,X,A,A', /silent
     wSN = where( obj->extract('clusterid') eq strmid( SNname, 3, 1 ) )
     if wSN[0] ne -1 then begin
        for iw=0, n_elements(wSN)-1 do begin
           if SNhostra[wSN[iw]] eq '??' then continue
           get_coords, coords, instring=SNhostra[wSN[iw]]+' '+SNhostdec[wSN[iw]]
           nearby = obj->nearestobjs(coords[0]*15., coords[1])
           print, SNname[wSN[iw]]+' = '+obj->extract('clusterid')+string(nearby[0], format='(I05)') + $
                  "mag = "+string(ra[nearby[0]], format='(F5.2)')+" color = "+string(colors[nearby[0]], format='(F5.2)')
           plotsym, 3, circsize*3, color=SNcolor
           oplot, [ra[nearby[0]]], [dec[nearby[0]]], ps=8
           xyouts, ra[nearby[0]], dec[nearby[0]], strmid( SNname[wSN[iw]], 3 ), $
                   charsize=2.0*circsize, charthick=2., align=-.5, color=SNcolor
        endfor
     endif
  endif
  

end
