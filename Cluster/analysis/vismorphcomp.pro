Pro vismorphcomp, sum, vm, color, hardcopy=hardcopy, title=title, plot2=plot2, morph=morph
  if n_elements(title) eq 0 then title=sum[0].clusterid+' morphology'

  if keyword_set(plot2) then begin
     xsize = 4.5
     ysize = 8.
     !p.multi=[0,1,2]
  endif else begin
     xsize = 4.5
     ysize = 4.5
  endelse

  if n_elements(hardcopy) ne 0 then begin
     thisDevice = !d.name
     set_plot, 'ps'
     device, /encapsulated, /color, bits_per_pixel=8, landscape=0
     device, xsize=xsize, ysize=ysize, /inches
     device, filename=hardcopy
     thick=4
     charthick=4
     xthick=4
     ythick=4
  endif else begin
     window, 0, xsize=100*xsize, ysize=100*ysize
     thick=1
     charthick=1
     xthick=1
     ythick=1
  endelse

  early = vm.ttype le -4
  s0 = vm.ttype ge -3 and vm.ttype le 0
  late = vm.ttype ge 1 and vm.ttype ne 99

  if n_elements(morph) eq 0 then morph=1
  case morph of 
     1 : morphcheck = sum.conc gt 0.35+4.65*sum.asym
     2 : morphcheck = sum.gini2 ge 0.44 and sum.asym2 le 0.08
  endcase

  red = color gt 0.85 and color lt 1.2

  wearlyfill = where(early and morphcheck and red)
  wearlyopen = where(early and ~morphcheck and red)

  ws0fill = where(s0 and morphcheck and red)
  ws0open = where(s0 and ~morphcheck and red)

  wlatefill = where(late and morphcheck and red)
  wlateopen = where(late and ~morphcheck and red)


  loadct, 3
  red = fsc_color('red', !d.table_size-2)
  purple = fsc_color('purple', !d.table_size-3)
  blue = fsc_color('blue', !d.table_size-4)

  plot, [0], [0], /nodata, xrange=[0.1, 0.7], yrange=[-0.05, 0.5], xstyle=1, ystyle=1, $
        xtitle='Gini2', ytitle='Asym2', charthick=charthick, $
        xthick=xthick, ythick=ythick, $
        title=title

  plotsym, 0, color=red, thick=thick, /fill
  if wearlyfill[0] ne -1 then $
     oplot, [sum[wearlyfill].gini2], [sum[wearlyfill].asym2], ps=8
  plotsym, 0, color=red, thick=thick
  if wearlyopen[0] ne -1 then $
     oplot, [sum[wearlyopen].gini2], [sum[wearlyopen].asym2], ps=8

  plotsym, 4, color=purple, thick=thick, /fill
  if ws0fill[0] ne -1 then $
     oplot, [sum[ws0fill].gini2], [sum[ws0fill].asym2], ps=8
  plotsym, 4, color=purple, thick=thick
  if ws0open[0] ne -1 then $
     oplot, [sum[ws0open].gini2], [sum[ws0open].asym2], ps=8

  plotsym, 8, color=blue, thick=thick, /fill
  if wlatefill[0] ne -1 then $
     oplot, [sum[wlatefill].gini2], [sum[wlatefill].asym2], ps=8
  plotsym, 8, color=blue, thick=thick
  if wlateopen[0] ne -1 then $
     oplot, [sum[wlateopen].gini2], [sum[wlateopen].asym2], ps=8


  if keyword_set(plot2) then begin
     plotsym, 0, color=red, thick=thick, /fill
     plot, [sum[wearlyfill].conc], [sum[wearlyfill].asym], ps=8, $
           xrange=[0.1, 0.8], yrange=[-0.05, 0.3], xstyle=1, ystyle=1, $
           xtitle='Conc', ytitle='Asym', charthick=charthick, $
           xthick=xthick, ythick=ythick
     plotsym, 0, color=red, thick=thick
     oplot, [sum[wearlyopen].conc], [sum[wearlyopen].asym], ps=8
     
     plotsym, 4, color=purple, thick=thick, /fill
     oplot, [sum[ws0fill].conc], [sum[ws0fill].asym], ps=8
     plotsym, 4, color=purple, thick=thick
     oplot, [sum[ws0open].conc], [sum[ws0open].asym], ps=8
     
     plotsym, 8, color=blue, thick=thick, /fill
     oplot, [sum[wlatefill].conc], [sum[wlatefill].asym], ps=8
     plotsym, 8, color=blue, thick=thick
     oplot, [sum[wlateopen].conc], [sum[wlateopen].asym], ps=8
     
  endif
  if n_elements(hardcopy) ne 0 then begin
     device, /close_file
     set_plot, thisDevice
  endif
  !p.multi=0
end
