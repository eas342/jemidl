Pro K_CMD, objs, _extra=extra
  resolve_routine,'clustergalaxies__define',/com, /no
  readcol, '/home/jmeyers314/supernova.txt', name, nickname, ra, dec, z, format='A,A,X,X,A,A,D,X,X,X,X,X'
  
  colors = Obj_New("IDLgrPalette")
  colors->LoadCT, 3
  colors->GetProperty, Red=r, Green=g, Blue=b
  Obj_Destroy, colors
  TVLCT, r, g, b
  grey = fsc_color('light gray', !d.table_size-4)
  red = fsc_color('red', !d.table_size-2)
  blue = fsc_color('blue', !d.table_size-3)
  orange = fsc_color('orange', !d.table_size-6)
  green = fsc_color('green', !d.table_size-7)


  thisDevice=!d.name
  set_plot, 'ps'
  device, /encapsulated, color=1, bits_per_pixel=8, landscape=0
  device, xsize=4.5, ysize=3., /inches
  device, filename='K_CMD.eps'
  
  obj=objs[9]
  colors = obj->color( _extra=extra )
  sexcat = obj->sexcat()
  speccat = obj->speccat()
  mags   = sexcat.zmag_best
  
  !p.font=1
  device, set_font='Times', /tt_font

  plot, [0], [0], /nodata, xstyle=1, ystyle=1, xrange=[19,25], yrange=[-1.,2.], $
        xtitle='z!D850', ytitle='i!D775!N - z!D850', title='ISCS J1438.1+3414 at z=1.41'

  ; K-cluster values...
  slope = -0.022740758
  int = 1.3429667
  disp=0.06393*3.0
  
  left=21
  right=25

  ul = (left)*slope+int + disp
  ll = (left)*slope+int - disp
  ur = (right)*slope+int + disp
  lr = (right)*slope+int - disp
  polyfill, [left, left, right, right, left], [ll, ul, ur, lr, ll], color=grey

  oplot, [left, right], [left, right]*slope+int

  w1 = where( speccat.z le 0. )
  plotsym, 0, 0.35
  oplot, mags[w1], colors[w1], ps=8

;  w2 = where( speccat.z ge 1.40 and speccat.z le 1.42 and colors lt 0.85 )
;  plotsym, 8, 0.45, /fill, color=red
;  oplot, mags[w2], colors[w2], ps=8

  w3 = where( speccat.z gt 0. and ~(speccat.z ge 1.40 and speccat.z le 1.42) )
  plotsym, 0, 0.35, /fill
  oplot, mags[w3], colors[w3], ps=8

;  w4 = where( speccat.z ge 1.40 and speccat.z le 1.42 )
;  plotsym, 0, 0.35, /fill
;  oplot, mags[w4], colors[w4], ps=8

  w = where( strmid( name, 3, 1 ) eq 'K' )
  for iw=0, n_elements(w)-1 do begin
     get_coords, coords, instring=ra[w[iw]]+' '+dec[w[iw]]
     alpha=coords[0]*15.d
     delta=coords[1]
     dist2 = (sexcat.alphawin_j2000 - alpha)^2*cos(delta*!dpi/180.d)^2 $
             + (sexcat.deltawin_j2000 - delta)^2
     md2 = min( dist2, m )*3600.d^2
     plotsym, 3, 0.85, /fill
     oplot, [mags[m]], [colors[m]], ps=8, color=green
  endfor

  xyouts, 21.1, 1.1, '?', charsize=1.1

  device, /close_file
  set_plot, thisDevice

end
