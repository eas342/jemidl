Pro measureRS, magrange=magrange, colorrange=colorrange, win=win1, wout=wout1, plotrange=plotrange, weightthresh=weightthresh, _extra=extra
  common JEM$cmdblock, colors, mags, $
     clustername, zcluster, clusterid, ra, dec, $
     size, ps, $
     wmorph, wnomorph, $
     wo2spec, wno2spec, wnocluespec, wnonspec, $
     wSN, SNname, $
     specmembers, specnonmembers, sne, $
     slope, intercept, sigma, $
     win, weight, wout, colorresid, z, $
     morphcolor

  if n_elements(weightthresh) eq 0 then weightthresh=0.8

  DoCMDPlot, _extra=extra
  red = fsc_color('red', !d.table_size-10)
  blue = fsc_color('blue', !d.table_size-11)
  
  if n_elements(win1) eq 0 then begin
     if n_elements(colorrange) eq 0 then begin
        case clusterid of 
           'A' : colorrange=[0.8, 1.0]
           'B' : colorrange=[0.9, 1.15]
           'C' : colorrange=[0.8, 1.05]
           'D' : colorrange=[0.85, 1.15]
           'E' : colorrange=[0.85, 1.15]
           'F' : colorrange=[0.9, 1.15]
           'G' : colorrange=[0.9, 1.15]
           'H' : colorrange=[0.85, 1.1]
           'I' : colorrange=[0.8, 1.1]
           'J' : colorrange=[0.8, 1.05]
           'K' : colorrange=[0.8, 1.05]
           'L' : colorrange=[0.8, 1.0]
           'M' : colorrange=[0.55, 0.95]
           'N' : colorrange=[0.8, 1.15]
           'P' : colorrange=[0.85, 1.1]
           'Q' : colorrange=[0.7, 1.05]
           'R' : colorrange=[0.85, 1.1]
           'S' : colorrange=[0.9, 1.1]
           'T' : colorrange=[0.75, 1.1]
           'U' : colorrange=[0.85, 1.1]
           'V' : colorrange=[0.6, 0.95]
           'W' : colorrange=[0.9, 1.15]
           'X' : colorrange=[0.85, 1.15]
           'Y' : colorrange=[0.85, 1.1]
           'Z' : colorrange=[0.75, 1.1]
           else : colorrange = [0.8,1.1]
        endcase

     endif
     if n_elements(magrange) eq 0 then begin
        case clusterid of 
           'A' : magrange = [22.7, 24]
           'B' : magrange = [20.6, 24]
           'C' : magrange = [20.5, 24]
           'D' : magrange = [20.6, 24]
           'E' : magrange = [20.5, 24]
           'F' : magrange = [21.2, 24]
           'G' : magrange = [22, 24]
           'H' : magrange = [21, 24]
           'I' : magrange = [21.5, 24]
           'J' : magrange = [21.5, 24]
           'K' : magrange = [22, 24]
           'L' : magrange = [22, 24]
           'M' : magrange = [20.2, 24]
           'N' : magrange = [20.8, 24]
           'P' : magrange = [21.5, 24]
           'Q' : magrange = [20.5, 24]
           'R' : magrange = [21, 24]
           'S' : magrange = [21.2, 24]
           'T' : magrange = [20.5, 24]
           'U' : magrange = [21, 24]
           'V' : magrange = [19.5, 24]
           'W' : magrange = [22, 24]
           'X' : magrange = [21.2, 24]
           'Y' : magrange = [21, 24]
           'Z' : magrange = [22, 24]
           else : magrange = [20, 24]
        endcase

     endif
     win = where(member(indgen(n_elements(mags)),wmorph) $
                 and colors ge colorrange[0] $
                 and colors le colorrange[1] $
                 and mags ge magrange[0] $
                 and mags le magrange[1])
     if keyword_set(plotrange) then begin
        oplot, [19,25], colorrange[0]*[1,1], color=blue
        oplot, [19,25], colorrange[1]*[1,1], color=blue
        oplot, magrange[0]*[1,1], [-0.5, 1.5], color=blue
        oplot, magrange[1]*[1,1], [-0.5, 1.5], color=blue
     endif
  endif else win=win1

  r=biweight_linfit(mags[win], colors[win], sigma, weight, /silent)
  slope=r[1]
  intercept=r[0]+r[1]*23.
  wh = where(weight gt weightthresh*max(weight))
  if wh[0] ne -1 then wout=win[wh]
  wout1=wout
  colorresid = colors-(r[0]+r[1]*mags)
;  plotsym, 0, size, /fill, color=red
;  oplot, mags[wout], colors[wout], ps=8
  oplot, [19,25],[19,25]*r[1]+r[0], color=red
  oplot, [19,25],[19,25]*r[1]+r[0]+sigma, linestyle=1, color=red
  oplot, [19,25],[19,25]*r[1]+r[0]-sigma, linestyle=1, color=red
  str = string(format='(%"b_23=%6.3f   m=%6.3f   s=%6.3f   n=%i")', intercept, slope, sigma, n_elements(wout))
  xyouts, 20.0, 1.2, str, charsize=size, charthick=1., align=0., color=morphcolor
  print, str
end

Pro Plot1CMDhist, _extra=extra
  common JEM$cmdblock
  w=where(member(indgen(n_elements(mags)),wmorph) and mags lt 24)
  c=colors[w]-slope*(mags[w]-25)
  plothist, c, xhist, yhist, peak=1, /noplot, xrange=[-0.5, 1.5], bin=0.05, /nan
;  plot, yhist, xhist, ps=10
  plotslopehist, xhist, yhist, slope, /yhist, yrange=[-0.5, 1.5], ystyle=1, xrange=[0,1.2], xstyle=1, color=morphcolor, xtickname=replicate(' ',30), ytickname=replicate(' ', 30), /noerase, _extra=extra
end

Pro RS_List, filename
  common JEM$cmdblock
  openw, lun, filename, /get_lun
  s=sort(mags[wout])
  
  for i=0, n_elements(wout)-1 do begin
     iobj=wout[s[i]]
     radec, ra[iobj], dec[iobj], ihr, imin, xsec, ideg, imn, xsc
     rastr = string(format='(%"%02i:%02i:%06.3f")', ihr, imin, xsec)
     decstr = string(format='(%"%+03i:%02i:%06.3f")', ideg, imn, xsc)
     printf, lun, format='(%"%6s %+7.3f %7.3f %s %s %6.4f")', clusterid+str(iobj), colorresid[iobj], $
             mags[iobj], rastr, decstr, z[iobj]
  endfor
  printf, lun
  w=where(abs(colorresid) le 0.1)
  s=sort(mags[w])
  for i=0, n_elements(w)-1 do begin
     iobj=w[s[i]]
     if member(iobj, wout) then continue
     radec, ra[iobj], dec[iobj], ihr, imin, xsec, ideg, imn, xsc
     rastr = string(format='(%"%02i:%02i:%06.3f")', ihr, imin, xsec)
     decstr = string(format='(%"%+03i:%02i:%06.3f")', ideg, imn, xsc)
     printf, lun, format='(%"%6s %+7.3f %7.3f %s %s %6.4f")', clusterid+str(iobj), colorresid[iobj], $
             mags[iobj], rastr, decstr, z[iobj]
  endfor

  close, lun
  free_lun, lun
end

Pro DoCMDPlot, _extra=extra
  common JEM$cmdblock

  tvlct, savered, savegreen, saveblue, /get
  if keyword_set(ps) then begin
     nomorphcolor = fsc_color('dark gray', !d.table_size-1)
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

  ;;;;;;;;;;;;;;;;;;;;
  ;; setup plot window
  
  plot, [0], [0], /nodata, xstyle=1, ystyle=1, $
        xrange=[19,25], yrange=[-0.5, 1.5], $
        title=string( format='(%"%s z=%5.2f")', $
                      clustername, $
                      zcluster), $
        color=morphcolor, $
        xtitle='z!D850!N (AB)', ytitle='i!D775!N-z!D850!N (AB)', $
        _extra=extra
  xyouts, 19.5, 1.1, clusterid, $
          charsize=4.0*size, charthick=2., align=0.5, color=morphcolor


  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ; plot all the gals...  will highlight next
  if wnomorph[0] ne -1 then begin
     plotsym, 0, size, color=nomorphcolor
     oplot, mags[wnomorph], colors[wnomorph], ps=8
  endif

  if wmorph[0] ne -1 then begin
     plotsym, 0, size, /fill, color=morphcolor
     oplot, mags[wmorph], colors[wmorph], ps=8
  endif

  if keyword_set(specmembers) then begin
     if wo2spec[0] ne -1 then begin
        plotsym, 4, size*2, color=o2color
        oplot, mags[wo2spec], colors[wo2spec], ps=8
     endif
     if wno2spec[0] ne -1 then begin
        plotsym, 8, size*2, color=no2color
        oplot, mags[wno2spec], colors[wno2spec], ps=8
     endif
     if wnocluespec[0] ne -1 then begin
        plotsym, 5, size*2, color=nocluecolor
        oplot, mags[wnocluespec], colors[wnocluespec], ps=8
     endif
  endif
  if keyword_set(specnonmembers) then begin
     if wnonspec[0] ne -1 then begin
        oplot, mags[wnonspec], colors[wnonspec], ps=7, symsize=size*1.8, color=noncolor
     endif
  endif
  if keyword_set(sne) then begin
     if n_elements(wSN) ne 0 then begin
        if wSN[0] ne -1 then begin
           for iw=0, n_elements(wSN)-1 do begin
              plotsym, 3, size*3, color=SNcolor
              oplot, [mags[wSN[iw]]], [colors[wSN[iw]]], ps=8
              xyouts, mags[wSN[iw]], -0.4+0.2*iw, strmid( SNname[iw], 3 ), $
                      charsize=2.0*size, charthick=2., align=0.5, color=SNcolor
           endfor
        endif
     endif
  endif

end

Pro Plot1CMD1ps, obj, filename=filename, _extra=extra
  thisDevice = !d.name
  set_plot, 'ps'
  device, /encapsulated, color=1, bits_per_pixel=8, /landscape
  device, xsize=10., ysize=8., /inches
  device, filename=filename
  plot1cmd1, obj, /ps, _extra=extra
  device, /close_file
  set_plot, thisDevice
end

Pro Plot1CMD1, obj, size=size1, ps=ps1, $
               specmembers=specmembers1, specnonmembers=specnonmembers1, sne=sne1, $
               hist=hist, morph=morph, _extra=extra
  common JEM$cmdblock

  slope=0
  intercept=1
  sigma=0.1
  resolve_obj, obj
  if n_elements(size1) eq 0 then size=1. else size=size1
  specmembers=keyword_set(specmembers1)
  specnonmembers=keyword_set(specnonmembers1)
  sne=keyword_set(sne1)
  ps=keyword_set(ps1)

  s=obj->summary()
  z=obj->extract('zcluster')
  colors = obj->color( _extra=extra )
  mags=s.zmag_auto
  z=s.z

  ebv = obj->extract('ebv')
  A_i = 1.973*ebv
  A_z = 1.472*ebv

  mags -= A_z
  colors -= (A_i - A_z)

  clustername = obj->extract('clustername')
  clusterid = obj->extract('clusterid')
  zcluster = obj->extract('zcluster')

  ;;;;;;;;;;;;;;;;;;;;;
  ; cut out things that are > 1 Mpc from center

  ra0=obj->extract('ra')
  dec0=obj->extract('dec')
  ra=s.alphawin_j2000
  dec=s.deltawin_j2000
  radius2 = (ra0-ra)^2*(cos(dec0*!dpi/180.d))^2 $
            + (dec0-dec)^2
  radius2 *= 3600.d^2
  radius_threshold = 1.0d/(lumdist( zcluster, /silent )/(1.+zcluster)^2)*180.d/!dpi*3600.d
  radcheck = radius2 lt radius_threshold^2

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ; cut out stars
  starcheck = s.zflux_radius lt 2.2

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ; cut out likely [oII] emitters
  if n_elements(morph) eq 0 then morph=1
  case morph of 
     1 : morphcheck = s.conc gt 0.35+4.65*s.asym
     2 : morphcheck = s.asym2 le 0.08 and s.gini2 gt 0.44
  endcase

;  morphcheck = s.conc gt 0.35+4.65*s.asym
  wmorph = where( morphcheck and radcheck and ~starcheck )
  wnomorph = where( ~(morphcheck and radcheck) and ~starcheck )
  ;;;;;;;;;;;;;;;;
  ; highlight spectroscopic members

  if keyword_set(specmembers) or keyword_set(specnonmembers) then begin
     
     if ~member('Z',tag_names(s)) then begin
        wno2spec = -1
        wo2spec = -1
        wnocluespec = -1
        wnonspec = -1
     endif else begin
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
           wno2spec = where(zcheck and no2check and (1-nocluecheck) and ~starcheck)
           wo2spec = where(zcheck and (1-no2check) and (1-nocluecheck) and ~starcheck)
           wnocluespec = where(zcheck and nocluecheck and ~starcheck)
        endif
        if keyword_set(specnonmembers) then begin
           wnonspec = where((1-zcheck) and s.z ge 0. and ~starcheck)
        endif
     endelse
  endif

  if keyword_set(sne) then begin
     readcol, '/home/scpdata02/clusters/supernova.txt', $
              SNename, SNnickname, SNhostra, SNhostdec, $
              format='A,A,X,X,X,X,X,A,A', /silent
     wSNe = where( obj->extract('clusterid') eq strmid( SNename, 3, 1 ) )
     delvarx, wSN, SNname
     if wSNe[0] ne -1 then begin
        for iw=0, n_elements(wSNe)-1 do begin
           if SNhostra[wSNe[iw]] eq '??' then continue
           get_coords, coords, instring=SNhostra[wSNe[iw]]+' '+SNhostdec[wSNe[iw]]
           nearby = obj->nearestobjs(coords[0]*15., coords[1])
           print, string(format='(%"%s = %s%05i  mags = %5.2f color = %5.2f asym2=%+7.3f gini2=%6.3f")', $
                         SNename[wSNe[iw]], obj->extract('clusterid'), nearby[0], $
                         mags[nearby[0]], colors[nearby[0]], s[nearby[0]].asym2, s[nearby[0]].gini2)
;           print, SNename[wSNe[iw]]+' = '+obj->extract('clusterid')+string(nearby[0], format='(I05)') + $
;                  "  mag = "+string(mags[nearby[0]], format='(F5.2)')+" color = "+string(colors[nearby[0]], format='(F5.2)')
           if n_elements(wSN) eq 0 then begin
              wSN = nearby[0]
              SNname = SNename[wSNe[iw]]
           endif else begin
              wSN=[wSN,nearby[0]]
              SNname = [SNname, SNename[wSNe[iw]]]
           endelse
        endfor
     endif
  endif
  if keyword_set(hist) then begin
     if total(!p.multi) gt 0 then begin
        region=pmultiposition()
        position=pmulti2position(region, [0.05,0.08,0.8,0.95])
        position2=pmulti2position(region, [0.8, 0.08, 0.95, 0.95])
        DoCMDPlot, position=position, /noerase, _extra=extra
        measurers, position=position, /noerase, _extra=extra
        plot1cmdhist, position=position2
     endif else begin
        position=[0.05,0.08,0.8,0.95]
        ;DoCMDPlot, position=position, _extra=extra
        measurers, position=position, _extra=extra
        position=[0.8, 0.08, 0.95, 0.95]
        plot1cmdhist, position=position
     endelse
  endif else DoCMDPlot, _extra=extra
end
