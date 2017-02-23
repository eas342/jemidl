Function MeasureO2MinimizeMe, $
   p

  common JEM$MeasureO2, bc300, flux, ivar, wcont, acoeff

  modelflux=bc300->flux( age=p[0], metal=p[1] )
  chi2 = computechi2( flux[wcont], $
                      sqrt(ivar[wcont]), $
                      modelflux[wcont], $
                      acoeff=acoeff )
  return, chi2
end

Function MeasureO2ModelFlux

  common JEM$MeasureO2

  mn = [2.5e9, 0.0177]
  for iter=0, 3 do begin ;;iterate
     mn = downhillsimplex( 'MeasureO2MinimizeMe', $
                           mn, $
                           [1.e8, 0.001] )
  endfor
  junk = measureo2minimizeme(mn) ;;update acoeff
  modelflux=acoeff[0]*bc300->flux( age=mn[0], metal=mn[1] )
  return, modelflux
end

Pro MeasureO2, $
   spec, $
   O2EWval=O2EWval, $
   O2EWerr=O2EWerr, $
   O2fluxval=O2fluxval, $
   O2fluxerr=O2fluxerr, $
   O2fluxul=O2fluxul, $
   plot=plot, $
   ulplot=ULplot
  

  O2EWval=!values.f_nan
  O2EWerr=!values.f_nan
  O2fluxval=!values.f_nan
  O2fluxerr=!values.f_nan
  O2fluxul=!values.f_nan

  common JEM$MeasureO2
  resolve_routine, 'spectrum__define', /no_recompile
  resolve_routine, 'dopplershift__define', /no_recompile

  if ~obj_valid(spec) then return
  
  ;; extract relavent data from spec object
  wave = spec->wavelength(frame='rest')
  dw = (wave-shift(wave,1))
  dw[0] = dw[1]
  flux = spec->flux(frame='rest')
  med = median(flux)
  flux /= med
;  ivar = spec->ivar(frame='rest')
  ivar = fltarr(n_elements(wave))+1.

  ;; define continuum bandpasses
  wblue = where( wave gt 3650. and wave lt 3700. )
  wred  = where( wave gt 3750. and wave lt 3800. )
  wcont = [wblue,wred]
  ;; O[II] bandpass
  wew = where( wave ge 3726.032*(1.-500./299792.) $
               and wave le 3728.815*(1.+500/299792.) )

  ;; punt if bandpasses arent possible
  if wblue[0] eq -1 or wred[0] eq -1 or wew[0] eq -1 then return
  if total(flux[wblue] ne 0.) eq 0 $
     or total(flux[wew] ne 0. ) eq 0 $
     or total(flux[wred] ne 0.) eq 0 then return
  if total(1-finite(flux[wblue])) ne 0 $
     or total(1-finite(flux[wew])) ne 0 $
     or total(1-finite(flux[wred])) ne 0 then return
  
  ;; create bc03 server blurred to 300km/s
  bc300 = obj_new( 'bc03server', wave=wave, sigmavbyc=0.001 )

  ;;now the actual measurement.
  modelflux=measureo2modelflux()
  if arg_present(O2EWerr) or arg_present(O2fluxerr) or arg_present(O2fluxul) then begin
     flux0 = flux
     wsig = biweight_scale( (modelflux-flux)[wcont], /zero )
     niter = 1000
     if arg_present(O2fluxul) then niter=1000
     ews = fltarr(niter)
     fluxs = fltarr(niter)
     for iter=0, niter-1 do begin
        counter, iter+1, niter, 'MC error analysis '
        flux = flux0 + randomn(seed, n_elements(flux))*wsig
        modelflux = measureo2modelflux()
        ews[iter] = total(((modelflux-flux)*dw/modelflux)[wew])
        fluxs[iter] = total(((flux-modelflux)*dw)[wew])
;        plot, wave, flux, xrange=[3600,3850], yrange=[-0.2, 0.2], xstyle=1, ystyle=1
;        wait, 0.01
     endfor
     print
     O2EWval = biweight_location( ews )
     O2EWerr = biweight_scale( ews-O2EWval, /zero )
     O2fluxval = biweight_location( fluxs )
     O2fluxerr = biweight_scale( fluxs-O2fluxval, /zero )
     O2fluxval *= med
     O2fluxerr *= med
     if arg_present(O2fluxul) then begin
        w=where(fluxs ge 0.)
        if w[0] eq -1 then o2fluxul = 0. $
        else o2fluxul = percentile( fluxs[w], 0.68 )*med
        if n_elements(ulplot) ne 0 then begin
           set_plot, 'ps'
           device, /encapsulated, /color, bits_per_pixel=8
           device, xsize=6, ysize=4, /inches, xoffset=0, yoffset=0
           device, filename=ulplot
           xrange = percentile( fluxs, [0.01,0.99] )
           diff = xrange[1]-xrange[0]
           xrange += [-1.,1.]*0.1*diff
           plothist, fluxs, xrange=xrange, yrange=[0,100], $
                     bin=(xrange[1]-xrange[0])/30, /nan, $
                     thick=3, xthick=2, ythick=2, charthick=2, $
                     xtitle='[OII] Luminosity ('+str(med)+')', $
                     ytitle='number', $
                     title=ulplot
           oplot, [0,0], [-1000,1000], thick=3
           oplot, [-100,100], [0,0], thick=3
           w=where(fluxs ge 0. and fluxs le o2fluxul/med)
           plothist, fluxs[w], $
                     bin=(xrange[1]-xrange[0])/30, /nan, $
                     thick=3, /fill, fcolor=fsc_color('green', 254), /fline, /over
           device, /close_file
           set_plot, 'X'
        endif
     endif
     flux=flux0
  endif else begin ;; simple case, no MC error analysis
     O2EWval = total(((modelflux-flux)*dw/modelflux)[wew])
     O2fluxval = total(((flux-modelflux)*dw)[wew])*med
  endelse

  if n_elements(plot) ne 0 then begin
     modelflux = measureo2modelflux()
     set_plot, 'ps'
     device, /encapsulated, /color, bits_per_pixel=8
     device, xsize=6, ysize=4, /inches, xoffset=0, yoffset=0
     device, filename=plot
     green=fsc_color('green', 254)
     red=fsc_color('red', 253)
     blue=fsc_color('blue', 252)
     plot, wave, flux*med, $
           xrange=[3600, 3850], thick=2, charthick=2, $
           xthick=2, ythick=2, xstyle=1, $
           xtitle='wavelength', $
           ytitle='Flux', $
           title=plot
     oplot, [3600, 3850], [0,0], thick=2
     if total(~finite(modelflux)) eq 0 then begin
        oplot, wave, modelflux*med, color=green, thick=2
        polyfill, [wave[wew],reverse(wave[wew])], $
                  [flux[wew], reverse(modelflux[wew])]*med, $
                  /line_fill, spacing=0.1, color=green, thick=2
        polyfill, [wave[wblue],reverse(wave[wblue])], $
                  [flux[wblue], reverse(modelflux[wblue])]*med, $
                  /line_fill, spacing=0.1, color=blue, thick=2
        polyfill, [wave[wred],reverse(wave[wred])], $
                  [flux[wred], reverse(modelflux[wred])]*med, $
                  /line_fill, spacing=0.1, color=red, thick=2
     endif
     o2str = string( format='(%"[OII] EW: %+9.2e +/- %8.2e")', $
                     O2EWval, O2EWerr )
     xyouts, 0.18, 0.85, o2str, /norm, charthick=2
     if n_elements(O2fluxval) ne 0 then begin
        o2fstr = string( format='(%"[OII] Flux: %+9.2e +/- %8.2e")', $
                         O2fluxval, O2fluxerr )
        xyouts, 0.18, 0.79, o2fstr, /norm, charthick=2
     endif
     device, /close_file
     set_plot, 'X'
  endif
  obj_destroy, bc300
  
end
