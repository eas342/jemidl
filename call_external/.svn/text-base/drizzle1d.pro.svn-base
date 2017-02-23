Pro drizzle1d, wave, $
               flux, $
               ivar=ivar, $
               outwave=outwave, $
               outflux=outflux, $
               outivar=outivar, $
               preflux=preflux, $
               preivar=preivar, $
               wavemin=wavemin, $
               wavemax=wavemax, $
               _extra=extra
; one of the following is required:
;              dlam=dlam, $
;              dloglam=dloglam, $
;              outbins=outbins


  if n_elements(wavemin) eq 0 then wavemin = min(wave)
  if n_elements(wavemax) eq 0 then wavemax = max(wave)

  ;; get outbins (won't overwrite if outbins is part of extra)
  outbins = build_outbins(wavemin=wavemin, wavemax=wavemax, _extra=extra)
  nbins = n_elements(outbins)

  ;; setup outflux/outivar
  nout = n_elements(outflux)
  if nout eq 0 then begin
     outflux = dblarr(nbins-1)
  endif else begin
     if nout ne nbins-1 then begin
        message, 'nbins-1 != nout!'
     endif
  endelse
  nivar = n_elements(outivar)
  if nivar eq 0 then begin
     outivar = dblarr(nbins-1)
  endif else begin
     if nivar ne nbins-1 then begin
        message, 'nbins-1 != nivar!'
     endif
  endelse
  ;; get input binning
  inbins = wave2bins(wave)

  ;; get inivar
  if n_elements(ivar) eq 0 then begin
     ivar = dblarr(n_elements(flux))+1.d
  endif

  if n_elements(preflux) eq 0 then preflux = dblarr(n_elements(outflux))
  if n_elements(preivar) eq 0 then preivar = dblarr(n_elements(outivar))

  ;; type casting
  inbins = double(inbins)
  flux = double(flux)
  ivar = double(ivar)
  outbins = double(outbins)
  preflux = double(preflux)
  preivar = double(preivar)

  ;;now we can call drizzle1d
  junk = call_external(get_jem_lib(), 'drizzle1d', $
                       flux, ivar, inbins, n_elements(inbins), $
                       outbins, n_elements(outbins), $
                       preflux, preivar, $
                       outflux, outivar, $
                       value=[0,0,0,1,0,1,0,0,0,0], /I_value, /CDECL, $
                       verbose=verbose, show_all_output=verbose)
  outwave = bins2wave(outbins)
end
