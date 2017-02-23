Function MeasureO2MinimizeMe, p
  common JEM$MeasureO2, bc300, flux, ivar, wcont, acoeff
  modelflux=bc300->flux( age=p[0], metal=p[1] )
  chi2 = computechi2( flux[wcont], $
                      sqrt(ivar[wcont]), $
                      modelflux[wcont], $
                      acoeff=acoeff )
  return, chi2
end

Function MeasureO2, spec, plot=plot, err=err
  common JEM$MeasureO2
  resolve_routine, 'spectrum__define', /no_recompile
  resolve_routine, 'dopplershift__define', /no_recompile

  if ~obj_valid(spec) then return, !values.f_nan
  

  wave = spec->wavelength(frame='rest')
  flux = spec->flux(frame='rest')*1.e17
;  ivar = spec->ivar(frame='rest')
  ivar = fltarr(n_elements(wave))+1.

  wcont = where( ( wave gt 3650. $
                   and wave lt 3700.) $
                 or ( wave gt 3750. $
                      and wave lt 3800. ) )
  if wcont[0] eq -1 then begin
     err=!values.f_nan
     return, !values.f_nan
  endif

  wew = where( wave ge 3726.032*(1-500./299792.) $
               and wave le 3728.815*(1+500/299792.) )
  if wew[0] eq -1 then begin
     err=!values.f_nan
     return, !values.f_nan
  endif
  
  bc300 = obj_new( 'bc03server', wave=wave, sigmavbyc=0.001 )
  mn = downhillsimplex( 'MeasureO2MinimizeMe', $
                        [2.5e9, 0.02], $
                        [1.e8, 0.001] )
  modelflux=acoeff[0]*bc300->flux( age=mn[0], metal=mn[1] )
  obj_destroy, bc300

  dw = (wave-shift(wave,1))[1:*]
  ew = total(((modelflux-flux)*dw/modelflux)[wew])
  if arg_present(err) then begin
     junk = biweight_mean((modelflux-flux)[wcont], wsig)
     wsig /= sqrt(n_elements(wew))
     ew1 = total(((modelflux+wsig-flux)*dw/modelflux)[wew])
     ew2 = total(((modelflux-wsig-flux)*dw/modelflux)[wew])
     ew3 = total(((modelflux-flux)*dw/(modelflux+wsig))[wew])
     ew4 = total(((modelflux-flux)*dw/(modelflux-wsig))[wew])
;     var=total((dw^2/ivar/modelflux^2)[wew])
;     err=sqrt(var)
     err=max(abs(ew-[ew1,ew2,ew3,ew4]))
  endif

  if n_elements(plot) ne 0 then begin
     set_plot, 'ps'
     device, /encapsulated, /color, bits_per_pixel=8
     device, xsize=6, ysize=4, /inches, xoffset=0, yoffset=0
     device, filename=plot
     green=fsc_color('green', 254)
     plot, wave, flux, xrange=[3650, 3800], thick=2, xthick=2, ythick=2, charthick=2, xstyle=1
     oplot, [3600, 3900], [0,0], thick=2
     if total(~finite(modelflux)) eq 0 then begin
        oplot, wave, modelflux, color=green, thick=2
        polyfill, [wave[wew],reverse(wave[wew])], $
                  [flux[wew], reverse(modelflux[wew])], $
                  /line_fill, spacing=0.1, color=green, thick=2
     endif
     o2str = string( format='(%"[OII] EW: %+8.2f +/- %6.2f")', ew, err )
     xyouts, 0.18, 0.85, o2str, /norm, charthick=2
     device, /close_file
     set_plot, 'X'
  endif

  return, ew

end
