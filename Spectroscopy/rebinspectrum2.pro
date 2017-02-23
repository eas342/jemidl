Function RebinSpectrum2, wave, flux, ivar=ivar, $
                         _extra=extra
; possible extras....
;                   dlam=dlam, $
;                   dloglam=dloglam, $
;                   outbins=outbins, $
;                   wavemin=wavemin, $
;                   wavemax=wavemax


  if n_elements(ivar) eq 0 then begin
     ivar = dblarr(n_elements(wave))+1.d
  endif
  drizzle1d, wave, flux, ivar=ivar, $
             outflux=outflux, outivar=outivar, $
             outwave=outwave, _extra=extra

  return, {wave:outwave, flux:outflux, ivar:outivar}
end
