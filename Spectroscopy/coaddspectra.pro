Function CoaddSpectra, wave1, flux1, wave2, flux2, $
                       ivar1=ivar1, ivar2=ivar2, $
                       wavemin=wavemin, wavemax=wavemax, $
                       _extra=extra
; necessary extras....
;                   dlam=dlam, $
;                   dloglam=dloglam

  if n_elements(ivar1) eq 0 then begin
     ivar1 = dblarr(n_elements(wave1))+1.d
  endif
  if n_elements(ivar2) eq 0 then begin
     ivar2 = dblarr(n_elements(wave2))+1.d
  endif

  if n_elements(outbins) eq 0 then begin
     if n_elements(wavemin) eq 0 then wavemin = min([wave1,wave2])
     if n_elements(wavemax) eq 0 then wavemax = max([wave1,wave2])
     outbins = build_outbins( wavemin=wavemin, wavemax=wavemax, _extra=extra )
  endif

  drizzle1d, wave1, flux1, ivar=ivar1, $
             outbins=outbins, $
             outflux=outflux, outivar=outivar
  preflux = outflux
  preivar = outivar
  drizzle1d, wave2, flux2, ivar=ivar2, $
             outbins=outbins, $
             preflux=preflux, preivar=preivar, $
             outwave=outwave, outflux=outflux, outivar=outivar
  return, {wave:outwave, flux:outflux, ivar:outivar}
end
