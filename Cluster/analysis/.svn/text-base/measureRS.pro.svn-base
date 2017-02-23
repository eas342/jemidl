Pro MeasureRS, obj, win=win, wout=wout, $
               colorlow=colorlow, colorhigh=colorhigh, $
               maglow=maglow, maghigh=maghigh, $
               _extra=extra

  red = fsc_color('red', !d.table_size-10)
  blue = fsc_color('blue', !d.table_size-11)
;  sexcat=obj->extract('sexcat')
  s=obj->summary()
  ebv = obj->extract('ebv')
  A_i = 1.973*ebv
  A_z = 1.472*ebv
  zmag = s.zmag_auto
  zmag -= A_z

  if n_elements(win) eq 0 then begin
     if n_elements(colorlow) eq 0 then colorlow = 0.8
     if n_elements(colorhigh) eq 0 then colorhigh = 1.1
     if n_elements(maglow) eq 0 then maglow = 20.
     if n_elements(maghigh) eq 0 then maghigh = 24.
     
     ra=obj->extract('ra')
     dec=obj->extract('dec')
     z=obj->extract('zcluster')
     radius2 = (ra - s.alphawin_j2000)^2*cos(dec*!dpi/180.d)^2 $
               + (dec - s.deltawin_j2000)^2
     radius2 *= 3600.d^2
     radius_threshold = 1.0d/(lumdist( z, /silent )/(1.+z)^2)*180.d/!dpi*3600.d
     radcheck = radius2 lt radius_threshold^2
     morphcat = obj->Extract('morphcat')
     morphcheck = morphcat.conc gt 0.36+4.5*morphcat.asym
     brightcheck = zmag gt maglow and zmag lt maghigh
     color = obj->color( _extra=extra )
     color -= (A_i-A_z)
     colorcheck = color gt colorlow and color lt colorhigh
     starcheck = s.zflux_radius lt 2.2
     w=where( morphcheck and radcheck and colorcheck and brightcheck and ~starcheck)
     color=color[w]
     zmag=zmag[w]
     if n_elements(win) ne 0 then w=win
     wout=w
  endif else begin
     color=obj->color( _extra=extra, gals=win )
     color -= (A_i-A_z)
     zmag=zmag[win]
  endelse
  w0=where(finite(color))
  r=biweight_linfit(zmag[w0],color[w0], sigma, weight)
  oplot, [19,25], colorlow*[1,1], color=blue
  oplot, [19,25], colorhigh*[1,1], color=blue
  oplot, maglow*[1,1], [-0.5, 1.5], color=blue
  oplot, maghigh*[1,1], [-0.5, 1.5], color=blue

  plotsym, 0, 1., /fill, color=red
  oplot, zmag, color, ps=8
  oplot, [19,25],[19,25]*r[1]+r[0], color=red
  oplot, [19,25],[19,25]*r[1]+r[0]+sigma, linestyle=1, color=red
  oplot, [19,25],[19,25]*r[1]+r[0]-sigma, linestyle=1, color=red
  

  print
  print, r[0]+23*r[1], r[1], sigma, total(weight/max(weight) gt 0.9)
end
