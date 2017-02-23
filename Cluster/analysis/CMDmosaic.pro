Pro CMDmosaic, objs, landscape=landscape, ps=ps, _extra=extra
  resolve_obj, objs[0]
  ps=keyword_set(ps)

  if keyword_set(ps) then begin
     thisDevice=!d.name
     set_plot, 'ps'
     device, /encapsulated, color=1, bits_per_pixel=8, landscape=0
     device, xsize=12., ysize=8., /inches
     if keyword_set(landscape) then begin
        device, /landscape
        device, xsize=10., ysize=8.0, /inches, xoffset=1, yoffset=1
     endif
     device, filename='CMD.eps'
  endif

  !p.multi=[0,5,5]

  zs=dblarr(n_elements(objs))
  for iobj=0, n_elements(objs)-1 do begin
     zs[iobj]=objs[iobj]->extract('zcluster')
  endfor
  sortz=sort(zs)

  for isort=0, n_elements(objs)-1 do begin
     z = zs[sortz[isort]]
     plot1CMD1, objs[sortz[isort]], size=0.25, ps=ps, _extra=extra
  endfor

  if keyword_set(ps) then begin
     device, /close_file
     set_plot, thisDevice
  endif
  !p.multi=0
  
end
