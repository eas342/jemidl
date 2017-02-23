Function DrizCoaddSpec, $
   spectra, $
   frame=frame, $
   rmserror=rmserror, $  ;; use rms error instead of statistical error 
   _extra=extra
  
  if n_elements(frame) eq 0 then frame='rest'

  wave=spectra[0]->wavelength(frame=frame)
  for ispec=1, n_elements(spectra)-1 do $
     wave = [wave,spectra[ispec]->wavelength(frame=frame)]
  wavemin=min(wave)
  wavemax=max(wave)

  outbins = build_outbins(wavemin=wavemin, wavemax=wavemax, _extra=extra)

  outflux = dblarr(n_elements(outbins)-1)
  outivar = dblarr(n_elements(outbins)-1)

  for ispec=0, n_elements(spectra)-1 do begin
     wave = spectra[ispec]->wavelength(frame=frame)
     flux = spectra[ispec]->flux(frame=frame, unit='flambda')
     ivar = spectra[ispec]->ivar(frame=frame, unit='flambda')
     if size(ivar, /tname) eq 'INT' then ivar=dblarr(n_elements(flux))+1.d

     preflux=outflux
     preivar=outivar
     jem_drizzle1d, wave, flux, ivar=ivar, $
                    outbins=outbins, $
                    preflux=preflux, preivar=preivar, $
                    outwave=outwave, outflux=outflux, outivar=outivar
  endfor

  if keyword_set(rmserror) then begin
     outfluxes = dblarr(n_elements(spectra), n_elements(outbins)-1)
     ncontrib = lonarr(n_elements(outbins)-1)
     for ispec=0, n_elements(spectra)-1 do begin
        wave = spectra[ispec]->wavelength(frame=frame)
        flux = spectra[ispec]->flux(frame=frame, unit='flambda')
        ivar = spectra[ispec]->ivar(frame=frame, unit='flambda')
        if size(ivar, /tname) eq 'INT' then ivar=dblarr(n_elements(flux))+1.d
        
        jem_drizzle1d, wave, flux, ivar=ivar, $
                       outbins=outbins, $
                       outflux=outflux1
        wzero = where(outflux1 eq 0, complement=wnonzero)
        if wzero[0] ne -1 then outflux1[wzero] = outflux[wzero]
        outfluxes[ispec,*] = (outflux1 - outflux)^2
        ncontrib[wnonzero] += 1
     endfor
     rms_sq = total(outfluxes, 1)/ncontrib
     outivar = 1./rms_sq
  endif

  outspec = obj_new('spectrum', outwave, outflux, ivar=outivar)
  return, outspec
end
