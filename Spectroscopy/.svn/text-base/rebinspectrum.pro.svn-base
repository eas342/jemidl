Function RebinSpectrum, spectrum, $
                        frame=frame, $
                        _extra=extra
; possible extras....
;                   dlam=dlam, $
;                   dloglam=dloglam, $
;                   outbins=outbins, $
;                   wavemin=wavemin, $
;                   wavemax=wavemax
  

  if n_elements(frame) eq 0 then frame = 'rest'
  
  wave = spectrum->wavelength(frame=frame)
  flux = spectrum->flux(frame=frame, unit='flambda')
  ivar = spectrum->ivar(frame=frame, unit='flambda')
  if size(ivar,/tname) eq 'INT' then begin
     noivar=1
     ivar = dblarr(n_elements(wave))+1.d
  endif
  jem_drizzle1d, wave, flux, ivar=ivar, $
                 outflux=outflux, outivar=outivar, $
                 outwave=outwave, _extra=extra
  
  new_z = obj_new('dopplershift', z=(spectrum->new_z())->z())
  if keyword_set(noivar) then begin
     outspec = obj_new('spectrum', outwave, outflux, $
                       unit='flambda', new_z=new_z )
  endif else begin
     outspec = obj_new('spectrum', outwave, outflux, $
                       ivar=outivar, unit='flambda', $
                       new_z=new_z) 
  endelse
  return, outspec
end
