Pro photcompare, obj, radius=radius
  loadct, 0
  blue=fsc_color('blue')
  resolve_obj, obj
  obj->Loadzimage
  sexcat=obj->extract('sexcat')
  galfitcat=obj->extract('galfitcat')
  zsex=sexcat.zmag_aper[9]
  sexback=sexcat.zbackground
  joshback=galfitcat.zsky
  !p.multi=[0,1,3]

  plot, zsex, joshback-sexback, ps=1, xrange=[19,26], yrange=[-0.002,0.002], xstyle=1, ystyle=1, xtitle='zmag', title='background difference'
  plot1stat=dblarr(6,2)
  for i=0,5 do begin
     w=where(zsex gt 20+i and zsex lt 21+i and abs(joshback-sexback) lt 0.005)
     mean = biweight_mean((joshback-sexback)[w],sigma)
     plot1stat[i,*] = [mean,sigma]
  endfor
  oplot, [19,26], [0,0]
  oplot, indgen(7)+20.5, plot1stat[*,0], color=blue
  plotsym, 8, color=blue
  oploterr, indgen(7)+20.5, plot1stat[*,0], plot1stat[*,1], 8
  xyouts, 20, 0.001, string(plot1stat[4,1])

  zphot=obj->apphot(5., band='z', /skysub)
  plot, zsex, zsex-zphot, ps=1, xrange=[19,26], yrange=[-0.1,0.1], xtitle='zmag', ytitle='joshmag vs sextractor mag'
  for i=0,5 do begin
     w=where(zsex gt 20+i and zsex lt 21+i and abs(zsex-zphot) lt 0.5 )
     mean = biweight_mean((zsex-zphot)[w],sigma)
     plot1stat[i,*] = [mean,sigma]
  endfor
  oplot, [19,26], [0,0]
  oplot, indgen(7)+20.5, plot1stat[*,0], color=blue
  plotsym, 8, color=blue
  oploterr, indgen(7)+20.5, plot1stat[*,0], plot1stat[*,1], 8
  xyouts, 20, 0.05, string(plot1stat[4,1])

  zphotsex=obj->apphot(5., band='z', /skysub, /sexsub)
  plot, zsex, zsex-zphotsex, ps=1, xrange=[19,26], yrange=[-0.1,0.1], xtitle='zmag', ytitle='joshmag with sextractor background vs sextractor mag'
  for i=0,5 do begin
     w=where(zsex gt 20+i and zsex lt 21+i and abs(zsex-zphotsex) lt 0.5)
     mean = biweight_mean((zsex-zphotsex)[w],sigma)
     plot1stat[i,*] = [mean,sigma]
  endfor
  oplot, [19,26], [0,0]
  oplot, indgen(7)+20.5, plot1stat[*,0], color=blue
  plotsym, 8, color=blue
  oploterr, indgen(7)+20.5, plot1stat[*,0], plot1stat[*,1], 8
  xyouts, 20, 0.05, string(plot1stat[4,1])

  !p.multi=0
end
