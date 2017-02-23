Pro oplot_rs, cmd, eps=eps, mega=mega, hist=hist, scale=scale, _extra=extra
  if n_elements(scale) eq 0 then scale=1.0
  junk = fsc_color(/allcolors, colorstructure=ctable)
  if keyword_set(eps) or keyword_set(mega) then begin
     thick=3
     color=ctable.black
  endif else begin ;;'X' display
     thick=1
     color=ctable.white
  endelse
  if ~ptr_valid(cmd.fit) then return
  rs = *cmd.fit
  
  oplot, [19,25], ([19,25]-23)*rs.slope+rs.intercept, color=ctable.red, thick=thick
  oplot, [19,25], ([19,25]-23)*rs.slope+rs.intercept+rs.mscatter, color=ctable.red, linestyle=1, thick=thick
  oplot, [19,25], ([19,25]-23)*rs.slope+rs.intercept-rs.mscatter, color=ctable.red, linestyle=1, thick=thick
  print, string(format='(%"b_23=%5.3f+/-%5.3f  m=%+6.3f+/-%5.3f  s=%5.3f+/-%5.3f  s_cor=%5.3f  n=%i")', $
                rs.intercept, rs.intercept_err, rs.slope, rs.slope_err, rs.mscatter, rs.mscatter_err, rs.iscatter, rs.ngals)
  if keyword_set(results) then begin
     xyouts, 20, 1.2, $
             string(format='(%"b_23=%5.3f+/-%5.3f  m=%+6.3f+/-%5.3f  s=%5.3f+/-%5.3f  s_cor=%5.3f  n=%i")', $
                    rs.intercept, rs.intercept_err, rs.slope, rs.slope_err, rs.mscatter, rs.mscatter_err, $
                    rs.iscatter, rs.n), charthick=thick
  endif

  if keyword_set(hist) then begin
     if total(!p.multi) gt 0 then begin
        !p.multi[0] += 1
        region=pmultiposition()
        trialregion = [0.9,0.95,0.9,0.92]+scale*[-0.10, -0.87, 0.05, 0.00]
        position=pmulti2position(region, trialregion)
     endif else begin
        trialregion = [0.9,0.95,0.9,0.92]+scale*[-0.10, -0.87, 0.05, 0.00]
        position=trialregion
     endelse

     xhist = *cmd.xhist
     yhist = *cmd.yhist
     plotslopehist, xhist, yhist, rs.slope*(.75/0.72), /yhist, $
                    yrange=[-0.5,1.5], ystyle=1, $
                    xrange=[0., 1.2], xstyle=1, $
                    color=color, $
                    xtickname=replicate(' ',30), ytickname=replicate(' ', 30), $
                    /noerase, position=position, $
                    xthick=thick, ythick=thick, thick=thick, _extra=extra
  endif

end

Pro plot_cmd, cmd, _extra=extra, eps=eps, $
              specmembers=specmembers, $
              specnonmembers=specnonmembers, $
              sne=sne, radlimit=radlimit, $
              rs=rs, hist=hist, plotregion=plotregion, $
              results=results, titlesize=titlesize, noid=noid, $
              noearly=noearly, nonoclue=nonoclue, xycharsize=xycharsize, $
              mega=mega, xycharthick=xycharthick, scale=scale, psize=psize

  if n_elements(scale) eq 0 then scale=1.0
  if n_elements(psize) eq 0 then psize=1.0
  
  junk=fsc_color(/allcolors, colorstructure=ctable)
  if n_elements(xycharsize) eq 0 then xycharsize=1.
  if keyword_set(eps) then begin
     nomorphcolor = ctable.slategray
     morphcolor = ctable.black
     background = ctable.white
     o2color = ctable.blue
     no2color = ctable.red
     noncolor = ctable.black
     SNcolor = ctable.orange
     nocluecolor = ctable.purple
     regioncolor = ctable.pink
     thick=3.
     charsize=3.*xycharsize
     charthick=5.
  endif else if keyword_set(mega) then begin
     nomorphcolor = ctable.slategray
     morphcolor = ctable.black
     background = ctable.white
     o2color = ctable.blue
     no2color = ctable.black
     noncolor = ctable.black
     SNcolor = ctable.red
     nocluecolor = ctable.purple
     regioncolor = ctable.pink
     thick=3.
     charsize=3.*xycharsize
     charthick=5.
  endif else begin ;;'X' display
     nomorphcolor = ctable.darkgray
     morphcolor = ctable.white
     background = ctable.black
     o2color = ctable.blue
     no2color = ctable.red
     noncolor = ctable.white
     SNcolor = ctable.gold
     nocluecolor = ctable.purple
     regioncolor = ctable.red8
     thick=1.
     charsize=2.*xycharsize
     charthick=1.
  endelse
  
  clustername = cmd.clustername
  clusterid   = cmd.clusterid
  zcluster    = cmd.zcluster

  if keyword_set(hist) then begin
     if total(!p.multi) gt 0 then begin
        if !p.multi[0] eq 0 then erase
        region=pmultiposition()
        trialregion=[0.9, 0.95, 0.9, 0.92]+scale*[-0.82, -0.87, -0.10, 0.00]
        position=pmulti2position(region, trialregion)
     endif else begin
        trialregion=[0.9, 0.95, 0.9, 0.92]+scale*[-0.82, -0.87, -0.10, 0.00]
        position=trialregion
        erase
     endelse
     extra=struct_addtags(extra, {position:position, $
                                  noerase:1})
  endif

  plot, [0], /nodata, xstyle=1, ystyle=1, $
        xrange=[19,25], yrange=[-0.5,1.5], $
        title=string( format='(%"%s z=%4.2f")', $
                      clustername, $
                      zcluster ), $
        xthick=thick, ythick=thick, $
        charthick=thick, $
        xtitle=textoidl('z_{850} (AB)'), $
        ytitle=textoidl('i_{775} - z_{850}'), $
        _extra=extra, charsize=titlesize
  if keyword_set(plotregion) then begin
     x = [cmd.magrange[0], cmd.magrange[0], cmd.magrange[1], cmd.magrange[1], cmd.magrange[0]]
     y = [cmd.colorrange[0], cmd.colorrange[1], cmd.colorrange[1], cmd.colorrange[0], cmd.colorrange[0]]
     y += cmd.slope*(x-23)
     polyfill, x, y, color = regioncolor
  endif

  if ~keyword_set(noid) then begin
     xyouts, 19.5, 1.1, clusterid, $
          charsize=charsize, charthick=thick, align=0.5, color=morphcolor
  endif

  gals=*cmd.gals

  if keyword_set(radlimit) then $
     wnomorph=where(~gals.morph and ~gals.star and gals.rad) $
  else $
     wnomorph=where(~gals.morph and ~gals.star or ~gals.rad)
  if wnomorph[0] ne -1 then begin
     plotsym, 0, psize, color=nomorphcolor, thick=thick
     oplot, [gals[wnomorph].zmag_auto], [gals[wnomorph].iz], ps=8
  endif
  wmorph=where(gals.morph and ~gals.star and gals.rad)
  if wmorph[0] ne -1 then begin
     plotsym, 0, psize, /fill, color=morphcolor, thick=thick
;     oplot, [gals[wmorph].zmag_auto], [gals[wmorph].iz], ps=8
     oploterror, [gals[wmorph].zmag_auto], [gals[wmorph].iz], [gals[wmorph].iz_err], ps=8, /nohat, errthick=thick
  endif

  if keyword_set(specmembers) then begin
     if keyword_set(radlimit) then $
        wo2spec = where(gals.o2 and gals.specmem and gals.rad) $
     else $
        wo2spec = where(gals.o2 and gals.specmem)
     if wo2spec[0] ne -1 and ~keyword_set(noearly) then begin
        plotsym, 4, 2.0*psize, color=o2color, thick=thick
        oplot, [gals[wo2spec].zmag_auto], [gals[wo2spec].iz], ps=8
     endif
     if keyword_set(radlimit) then $
        wno2spec = where(gals.no2 and gals.specmem and gals.rad) $
     else $
        wno2spec = where(gals.no2 and gals.specmem)        
     if wno2spec[0] ne -1 then begin
        plotsym, 8, 2.0*psize, color=no2color, thick=thick
        oplot, [gals[wno2spec].zmag_auto], [gals[wno2spec].iz], ps=8
     endif
     if keyword_set(radlimit) then $
        wnocluespec = where(gals.noclue and gals.specmem and gals.rad) $
     else $
        wnocluespec = where(gals.noclue and gals.specmem)
     if wnocluespec[0] ne -1 and ~keyword_set(nonoclue) then begin
        plotsym, 5, 2.0*psize, color=nocluecolor, thick=thick
        oplot, [gals[wnocluespec].zmag_auto], [gals[wnocluespec].iz], ps=8
     endif
  endif
  
  if keyword_set(specnonmembers) then begin
     if keyword_set(radlimit) then $
        wnonspec = where(gals.specnon and gals.rad) $
     else $
        wnonspec = where(gals.specnon)
     if wnonspec[0] ne -1 then begin
        oplot, [gals[wnonspec].zmag_auto], [gals[wnonspec].iz], ps=7, symsize=2.0*psize, color=noncolor, thick=thick
     endif
  endif

  if keyword_set(sne) then begin
     w=where(gals.snname ne '' and gals.sntype eq 'SNIa')
     if w[0] ne -1 then begin
        for iw=0, n_elements(w)-1 do begin
           print, string(format='(%"%s = %s%05i  zmag = %5.2f color = %5.2f")', $
                         gals[w[iw]].snname, clusterid, gals[w[iw]].galid, $
                         gals[w[iw]].zmag_auto, gals[w[iw]].iz )
           plotsym, 3, 3.0*psize, color=SNcolor, thick=thick
           oplot, [gals[w[iw]].zmag_auto], [gals[w[iw]].iz], ps=8
           xyouts, gals[w[iw]].zmag_auto, -0.4+0.2*iw, gals[w[iw]].snname, $
                   charsize=charsize, charthick=charthick, align=0.5, color=SNcolor
        endfor
     endif
  endif
  if keyword_set(results) then begin
     fit=*cmd.fit
     xyouts, 19.2, 1.3, $
             string(format='(%"b_23=%5.3f+/-%5.3f  m=%+6.3f+/-%5.3f  s=%5.3f+/-%5.3f  s_cor=%5.3f  n=%i")', $
                    fit.intercept, fit.intercept_err, fit.slope, fit.slope_err, fit.mscatter, $
                    fit.mscatter_err, fit.iscatter, fit.ngals), $
             charsize=1., charthick=2
  endif
  if keyword_set(rs) then begin
     oplot_rs, cmd, hist=keyword_set(hist), eps=keyword_set(eps), mega=keyword_set(mega), scale=scale
  endif
end
