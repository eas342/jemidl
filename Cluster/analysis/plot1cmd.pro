Pro Plot1CMD, obj, circsize, ps=ps, specmembers=specmembers, specnonmembers=specnonmembers, sne=sne, _extra=extra
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
  endif else begin
     nomorphcolor = fsc_color('dark gray', !d.table_size-1)
     morphcolor = fsc_color('white', !d.table_size-2)
     background = fsc_color('black', !d.table_size-3)
     o2color = fsc_color('blue', !d.table_size-4)
     no2color = fsc_color('red', !d.table_size-5)
     noncolor = fsc_color('white', !d.table_size-6)
     SNcolor = fsc_color('gold', !d.table_size-7)
     nocluecolor = fsc_color('purple', !d.table_size-8)
  endelse

  s=obj->summary()
  z=obj->extract('zcluster')
  colors = obj->color( _extra=extra )
  mags = s.zmag_auto

  ebv = obj->extract('ebv')
  A_i = 1.973*ebv
  A_z = 1.472*ebv
  
  mags -= A_z

  colors -= (A_i - A_z)

  ;;;;;;;;;;;;;;;;;;;; 
  ; setup plot window

  plot, [0], [0], /nodata, xstyle=1, ystyle=1, $
        xrange=[19,25], yrange=[-0.5, 1.5], $
        title=string( format='(%"%s z=%5.2f")', obj->extract('clustername'), obj->extract('zcluster')), $
        color=morphcolor, xtitle='z!D850!N (AB)', ytitle='i!D775!N-z!D850!N (AB)'
  xyouts, 19.5, 1.1, obj->extract('clusterid'), charsize=4.0*circsize, charthick=2., align=0.5
  ra=obj->extract('ra')
  dec=obj->extract('dec')
  
  ;;;;;;;;;;;;;;;;;;;;;
  ; cut out things that are > 1 Mpc from center

  radius2 = (ra-s.alphawin_j2000)^2*(cos(dec*!dpi/180.d))^2 $
            + (dec-s.deltawin_j2000)^2
  radius2 *= 3600.d^2
  radius_threshold = 1.0d/(lumdist( z, /silent )/(1.+z)^2)*180.d/!dpi*3600.d
  radcheck = radius2 lt radius_threshold^2


  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ; cut out stars
  starcheck = s.zflux_radius lt 2.2

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ; cut out likely [oII] emitters

  morphcheck = s.conc gt 0.36+4.5*s.asym
  wmorph = where( morphcheck and radcheck and ~starcheck )
  wnomorph = where( ~(morphcheck and radcheck) and ~starcheck )
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ; plot all the gals...  will highlight next

  if wnomorph[0] ne -1 then begin
     plotsym, 0, circsize, color=nomorphcolor
     oplot, mags[wnomorph], colors[wnomorph], ps=8
  endif

  if wmorph[0] ne -1 then begin
     plotsym, 0, circsize, /fill, color=morphcolor
     oplot, mags[wmorph], colors[wmorph], ps=8
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
        'U' : zcheck = s.z ge 1.02  and s.z le 1.045
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
        wno2spec = where(zcheck and no2check and (1-nocluecheck) and ~starcheck)
        wo2spec = where(zcheck and (1-no2check) and ~starcheck)
        wnocluespec = where(zcheck and nocluecheck and ~starcheck)
        if wo2spec[0] ne -1 then begin
           plotsym, 8, circsize*2, color=o2color
           oplot, mags[wo2spec], colors[wo2spec], ps=8
        endif
        if wno2spec[0] ne -1 then begin
           plotsym, 8, circsize*2, color=no2color
           oplot, mags[wno2spec], colors[wno2spec], ps=8
        endif
        if wnocluespec[0] ne -1 then begin
           plotsym, 8, circsize*2, color=nocluecolor
           oplot, mags[wnocluespec], colors[wnocluespec], ps=8
        endif
     endif
     if keyword_set(specnonmembers) then begin
        wnonspec = where((1-zcheck) and s.z ge 0. and ~starcheck)
        if wnonspec[0] ne -1 then begin
           oplot, mags[wnonspec], colors[wnonspec], ps=7, symsize=circsize*1.8, color=noncolor
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
                  "  mag = "+string(mags[nearby[0]], format='(F5.2)')+" color = "+string(colors[nearby[0]], format='(F5.2)')
           plotsym, 3, circsize*3, color=SNcolor
           oplot, [mags[nearby[0]]], [colors[nearby[0]]], ps=8
           xyouts, mags[nearby[0]], -0.4+0.2*iw, strmid( SNname[wSN[iw]], 3 ), $
                   charsize=2.0*circsize, charthick=2., align=0.5, color=SNcolor
        endfor
     endif
  endif
end
