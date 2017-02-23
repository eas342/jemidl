Pro testdrizzle
  readcol, getenv('JEM_REF')+'SPECTRA/vega.dat', wave, flux
  drizzle1d, wave, flux, dloglam=0.003, $
             outwave=outwave, outflux=outflux
  splot, wave, flux
  soplot, outwave, outflux, color='blue'
  print, n_elements(wave)
  print, n_elements(outwave)
end

Pro testdrizzle2
  flux=[1.0d]
  ivar=[1.0d]
  wave=[1.0d]
  inbins=[0.0d, 1.0d]
  outbins=[0.0d, 0.25d, 0.75d, 1.0d]
  outflux = dblarr(n_elements(outbins)-1)
  outivar = dblarr(n_elements(outbins)-1)
  preflux = dblarr(n_elements(outflux))
  preivar = dblarr(n_elements(outivar))
  junk = call_external(get_jem_lib(), 'drizzle1d', $
                       flux, ivar, inbins, n_elements(inbins), $
                       outbins, n_elements(outbins), $
                       preflux, preivar, $
                       outflux, outivar, $
                       value=[0,0,0,1,0,1,0,0,0,0], /I_value, /CDECL, $
                       verbose=verbose, show_all_output=verbose)
  outflux2 = dblarr(1)
  outivar2 = dblarr(1)
  junk = call_external(get_jem_lib(), 'drizzle1d', $
                       outflux, outivar, outbins, n_elements(outbins), $
                       inbins, n_elements(inbins), $
                       preflux, preivar, $
                       outflux2, outivar2, $
                       value=[0,0,0,1,0,1,0,0,0,0], /I_value, /CDECL, $
                       verbose=verbose, show_all_output=verbose)
  stop
end
