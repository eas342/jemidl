Pro colors, background=background, _extra=extra
  junk=fsc_color(/all, color=ctable)
  plot, [0], /nodata, xrange=[0,1], yrange=[0,1], color=ctable.black

  if n_elements(background) ne 0 then begin
     w=where(strmatch(tag_names(ctable), background, /fold_case))
     polyfill, [0,1,1,0,0], [0,0,1,1,0], $
               color=ctable.(w), /normal
  endif
  tags = tag_names(ctable)
  for i=0, n_elements(tags)-1 do begin
     x=i / 40
     y=i mod 40
     xyouts, x/5., y/40., tags[i], color=ctable.(i), _extra=extra
  endfor
end
