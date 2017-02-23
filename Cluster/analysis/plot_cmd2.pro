Pro oplot_rs2, cmd, eps=eps, mega=mega, hist=hist, scale=scale, results=results, _extra=extra
  if n_elements(scale) eq 0 then scale=1.0
  junk = fsc_color( /all, color=ctable )
  if keyword_set(eps) or keyword_set(mega) then begin
     thick=3
     color=ctable.black
  endif else begin
     thick=1
     color=ctable.white
  endelse
  if ~ptr_valid(cmd.fit) then return
  rs = *cmd.fit

  if ptr_valid((*cmd.fit).xplot) then begin
     oplot, [19,25], ([19,25]-23)*rs.slope+rs.intercept, $
            color=ctable.red, thick=thick
     oplot, [19,25], ([19,25]-23)*rs.slope+rs.intercept+rs.mscatter, $
            color=ctable.red, linestyle=1, thick=thick
     oplot, [19,25], ([19,25]-23)*rs.slope+rs.intercept-rs.mscatter, $
            color=ctable.red, linestyle=1, thick=thick

     print, string( format='(%"b23=%5.3f  m=%+6.3f  mscat=%5.3f  cor=%5.3f  iscat=%5.3f")', $
                    rs.intercept, rs.slope, rs.mscatter, rs.scatter_corr, abs(rs.iscatter) )
  endif

  if keyword_set(results) then begin
     xyouts, 19.5, 1.35, $
             string( format='(%"b23=%5.3f  m=%+6.3f  mscat=%5.3f' $
                     +'  cor=%5.3f  iscat=%5.3f  n=%4.1f")', $
                     rs.intercept, rs.slope, rs.mscatter, $
                     rs.scatter_corr, abs(rs.iscatter), rs.ngals ), charthick=thick
     
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

     xhist = *(*cmd.fit).xhist
     yhist = *(*cmd.fit).yhist
     center = (*cmd.fit).intercept+(*cmd.fit).slope*(25-23)
     plot, [0], [0], xrange=[0., 1.2], yrange=[-0.5,1.5], xstyle=1, ystyle=1, $
           xtickname=replicate(' ',30), ytickname=replicate(' ', 30), $
           /noerase, position=position, $
           xthick=thick, ythick=thick, thick=thick, _extra=extra, color=ctable.black
           
     bin = xhist[1]-xhist[0]
     plotslopehist, xhist+center-bin/2., yhist, rs.slope*(0.75/0.72), /yhist, $
                    yrange=[-0.5,1.5], ystyle=5, $
                    xrange=[0., 1.2], xstyle=5, $
                    color=ctable.red, $
                    xtickname=replicate(' ',30), ytickname=replicate(' ', 30), $
                    /noerase, position=position, $
                    xthick=thick, ythick=thick, thick=thick, _extra=extra
     if ptr_valid((*cmd.fit).bkgxhist) then begin
        xhist = *(*cmd.fit).bkgxhist
        yhist = *(*cmd.fit).bkgyhist
        bin = xhist[1]-xhist[0]
        plotslopehist, xhist+center-bin/2., yhist, rs.slope*(0.75/0.72), /yhist, $
                       color=ctable.blue, $
                       xrange=[0.,1.2], yrange=[-0.5,1.5], xstyle=5, ystyle=5, $
                       xtickname=replicate(' ',30), ytickname=replicate(' ', 30), $
                       /noerase, position=position, $
                       xthick=thick, ythick=thick, thick=thick, _extra=extra
     endif
     if ptr_valid((*cmd.fit).xplot) and ~member(cmd.clusterid,['J','L']) then begin
        yplot = *(*cmd.fit).xplot + rs.slope*(25-23)
        xplot = *(*cmd.fit).yplot/(*cmd.fit).histmax
        yplot += rs.slope*(0.75/0.72)*xplot
        oplot, xplot, yplot, $
               thick=thick, color=ctable.magenta
     endif
     if ptr_valid((*cmd.fit).bkgxplot) and ~member(cmd.clusterid,['J','L']) then begin
        yplot = *(*cmd.fit).bkgxplot + rs.slope*(25-23)
        xplot = *(*cmd.fit).bkgyplot/(*cmd.fit).histmax
        yplot += rs.slope*(0.75/0.72)*xplot
        oplot, xplot, yplot, $
               thick=thick, color=ctable.cyan
     endif
  endif


end

Pro plot_cmd2, cmd, eps=eps, specmembers=specmembers, $
               specnonmembers=specnonmembers, $
               sne=sne, radlimit=radlimit, rs=rs, $
               hist=hist, plotselect=plotselect, $
               results=results, titlesize=titlesize, $
               noid=noid, nolate=nolate, nonoclue=nonoclue, $
               xycharsize=xycharsize, mega=mega, $
               xycharthick=xycharthick, scale=scale, psize=psize, $
               ext_arrow=ext_arrow, _extra=extra

  if n_elements(scale) eq 0 then scale=1.0
  if n_elements(psize) eq 0 then psize=1.0
  if n_elements(xycharsize) eq 0 then xycharsize=1.
  if n_elements(titlesize) eq 0 then titlesize=1.

  junk=fsc_color(/allcolors, colorstructure=ctable)
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

  ;; make room for histogram if requested
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

  ;;make initial blank plot
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

  ;;write out cluster id if not disabled
  if ~keyword_set(noid) then begin
     xyouts, 19.5, 1.1, clusterid, $
             charsize=charsize, charthick=thick, align=0.5, color=morphcolor
  endif

  gals=*cmd.gals
  
  ;;plot unselected galaxies
  wnoselect = where(~gals.select and ~gals.starcheck and gals.radcheck)
  plotsym, 0, psize, color=nomorphcolor, thick=thick
  oplot, [gals[wnoselect].z850], [gals[wnoselect].iz], ps=8
  ;;plot selected galaxies
  wselect = where(gals.select and ~gals.starcheck)
  plotsym, 0, psize, /fill, color=morphcolor, thick=thick
  oploterror, [gals[wselect].z850], $
              [gals[wselect].iz], $
              [gals[wselect].izerr], $
              ps=8, /nohat, errthick=thick
  
  ;;highlight spectroscopic members
  if keyword_set(specmembers) then begin
     wo2spec = where(gals.o2 and gals.specmem)
     if wo2spec[0] ne -1 and ~keyword_set(nolate) then begin
        plotsym, 4, 2.0*psize, color=o2color, thick=thick
        oplot, [gals[wo2spec].z850], [gals[wo2spec].iz], ps=8
     endif
     wno2spec = where(gals.no2 and gals.specmem)
     if wno2spec[0] ne -1 then begin
        plotsym, 8, 2.0*psize, color=no2color, thick=thick
        oplot, [gals[wno2spec].z850], [gals[wno2spec].iz], ps=8
     endif
     wnocluespec = where(gals.noclue and gals.specmem)
     if wnocluespec[0] ne -1 and ~keyword_set(nonoclue) then begin
        plotsym, 5, 2.0*psize, color=nocluecolor, thick=thick
        oplot, [gals[wnocluespec].z850], [gals[wnocluespec].iz], ps=8
     endif
  endif
  ;;highlight spectroscopic nonmembers
  if keyword_set(specnonmembers) then begin
     wnonspec = where(gals.specnon)
     if wnonspec[0] ne -1 then begin
        oplot, [gals[wnonspec].z850], [gals[wnonspec].iz], $
               ps=7, symsize=2.0*psize, color=noncolor, thick=thick
     endif
  endif

  ;;highlight SNe
  if keyword_set(sne) then begin
     w=where(gals.snname ne '' and gals.sntype eq 'SNIa')
     if w[0] ne -1 then begin
        for iw=0, n_elements(w)-1 do begin
           print, string(format='(%"%s = %s%05i  zmag = %5.2f color = %5.2f")', $
                         gals[w[iw]].snname, clusterid, gals[w[iw]].galid, $
                         gals[w[iw]].z850, gals[w[iw]].iz )
           plotsym, 3, 3.0*psize, color=SNcolor, thick=thick
           oplot, [gals[w[iw]].z850], [gals[w[iw]].iz], ps=8
           xyouts, gals[w[iw]].z850, -0.4+0.2*iw, gals[w[iw]].snname, $
                   charsize=charsize, charthick=charthick, align=0.5, color=SNcolor
        endfor
     endif
  endif

  if keyword_set(plotselect) and ptr_valid(cmd.fit) then begin
     if ptr_valid((*cmd.fit).inliers) then begin
        inliers = *(*cmd.fit).inliers
        plotsym, 0, psize*2.0, color=morphcolor, thick=thick
        oplot, [gals[inliers].z850], [gals[inliers].iz], ps=8
     endif
     if ptr_valid((*cmd.fit).outliers) then begin
        outliers = *(*cmd.fit).outliers
        plotsym, 0, psize*2.0, color=ctable.red, thick=thick
        oplot, [gals[outliers].z850], [gals[outliers].iz], ps=8
     endif
     oplot, [1.,1.]*(*cmd.fit).magrange[0], [-0.5, 0.0], color=morphcolor, thick=thick
     oplot, [1.,1.]*(*cmd.fit).magrange[1], [-0.5, 0.0], color=morphcolor, thick=thick
  endif

  if keyword_set(ext_arrow) then begin
     spec = bc03spec( 2.5e9, 0.017, /erg, origmass=1.e11 )
     wave = spec->wavelength()
     lumdens = spec->flux()
     iz0 = izmags4( wave, lumdens, cmd.zcluster, host_ebv=0. )
     iz1 = izmags4( wave, lumdens, cmd.zcluster, host_ebv=1./3.1 )
     color0 = iz0[0]-iz0[1]
     color1 = iz1[0]-iz1[1]
     dcolor = color1-color0
     dmag = iz1[1]-iz0[1]
     arrow, 20, -0.3, 20+dmag, -0.3+dcolor, thick=thick, color=morphcolor, /data
     xyouts, 20.0, -0.15, textoidl('A_V = 1.0'), charthick=2
  endif

  if keyword_set(rs) then begin
     oplot_rs2, cmd, hist=keyword_set(hist), eps=keyword_set(eps), $
                mega=keyword_set(mega), results=keyword_set(results), scale=scale
  endif
end
