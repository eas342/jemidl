Pro sensfuncplot, stds, sens, calspec, $
                  msk_balm=msk_balm1
  msk_balm=keyword_set(msk_balm1)
  nstd=n_elements(stds)
  sensfitarr = ptrarr(nstd)
  wavearr = ptrarr(nstd)
  for i=0, nstd-1 do begin
     std = xmrdfits(stds[i], 5, /silent)
     nspec = n_elements(stds)
     sig=fltarr(nspec)
     for j=0, nspec-1 do begin
        sig[j] = median(std[j].flux_opt)
     endfor
     junk = max(sig, m)
     std=std[m]
     sens1 = long_sensfunc(stds[i], sens[i], std_name=calspec[i], $
                           msk_balm=msk_balm, $
                           wave=wave, sensfunc=sensfunc, sensfit=sensfit)
     sensfitarr[i] = ptr_new(sensfit)
     wavearr[i] = ptr_new(wave)

     longslit_dir = getenv('LONGSLIT_DIR')
     cal_file = longslit_dir+'/calib/standards/calspec/'+calspec[i]+'.fits.gz'
     cal = xmrdfits(cal_file, 1, /silent)
     long_fluxcal, stds[i], /frm_sci, sensfuncfile=sens[i], wave=wavefx, flux=fluxfx

     ps_start, filename=repstr(sens[i],'.fits', '.eps'), /encapsulated, /nomatch
     xrange = minmax(std.wave_opt)
     junk = fsc_color(/all, color=ctable, /check_connection)
     multiplot, [1,4], mTitle=sens[i]+' '+calspec[i]
     yrange = percentile(std.flux_opt, [0.01, 0.99])*[1./1.1, 1.1]
     plot, /ylog, std.wave_opt, std.flux_opt, $
           xrange=xrange, xstyle=1, yrange=yrange, /ystyle

     multiplot
     w=where(cal.wavelength gt xrange[0] and cal.wavelength lt xrange[1])
     yrange = percentile(cal[w].flux, [0.01, 0.99])*[1./1.1, 1.1]
     plot, /ylog, cal.wavelength, cal.flux, $
           xrange=xrange, xstyle=1, yrange=yrange, /ystyle

     multiplot
     yrange = percentile(sensfit, [0.01, 0.99])*[1./1.1, 1.1]
     plot, /ylog, wave, sensfunc, $
           xrange=xrange, xstyle=1, yrange=yrange, /ystyle
     oplot, wave, sensfit, color=ctable.red

     multiplot, /doxaxis
     w=where(cal.wavelength gt xrange[0] and cal.wavelength lt xrange[1])
     yrange = percentile(cal[w].flux, [0.01, 0.99])*[1./1.1, 1.1]
     plot, /ylog, cal.wavelength, cal.flux, $
           xrange=xrange, xstyle=1, yrange=yrange, /ystyle
     oplot, wavefx, fluxfx/1.d17, color=ctable.red

     multiplot, /reset
     ps_end
  endfor
  xrange = minmax(*wavearr[0])
  yrange = percentile(*sensfitarr[0], [0.01, 0.99])
  ;; normalize to median wavelength
  junk=min(abs(*wavearr[0] - median(*wavearr[0])), w)
  wave0 = (*wavearr[0])[w]
  sens0 = (*sensfitarr[0])[w]
  for i=0, nstd-1 do begin
     if xrange[0] gt min(*wavearr[i]) then xrange[0] = min(*wavearr[i])
     if xrange[1] lt max(*wavearr[i]) then xrange[1] = max(*wavearr[i])
     sens1 = interpol((*sensfitarr[i]), (*wavearr[i]), wave0)
     p = percentile(*sensfitarr[i]/sens1*sens0, [0.01, 0.99])
     if yrange[0] gt p[0] then yrange[0] = p[0]
     if yrange[1] lt p[1] then yrange[1] = p[1]
  endfor
  yrange *= [1./1.1, 1.1]

  colors=[ctable.black, ctable.red, ctable.blue, ctable.magenta, ctable.cyan, ctable.green, ctable.yellow, ctable.orange]
  ps_start, filename='sensfunc_comparison.eps', /encapsulated, /nomatch
  plot, [0], /nodata, xrange=xrange, yrange=yrange, /xstyle, /ystyle, /ylog
  for i=0, nstd-1 do begin
     sens1 = interpol((*sensfitarr[i]), (*wavearr[i]), wave0)
     oplot, *wavearr[i], *sensfitarr[i]/sens1*sens0, color=colors[i]
  endfor
;  al_legend, sens+' '+calspec, colors=colors[0:nstd-1], textcolors=colors[0:nstd-1]
  ps_end
end
