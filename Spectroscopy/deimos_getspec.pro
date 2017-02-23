pro deimos_getspec, filename, wave, flux, ivar, wave_air=wave_air, $
                    sensfunc=sensfunc, senswave=senswave, $
                    raw=raw, opt=opt, box=box
  if (n_elements(opt) eq 0) then opt=0
  if (opt) then begin
     b = mrdfits(filename, 3, /silent)
     r = mrdfits(filename, 4, /silent)
  endif else begin
     b = mrdfits(filename, 1, /silent)
     r = mrdfits(filename, 2, /silent)
  endelse
  if n_elements(raw) eq 0 then raw=0
  if(raw) then begin
     wave = [b.lambda, r.lambda]
     flux = [b.spec, r.spec]
     ivar = [b.ivar, r.ivar]
  endif else begin
     wave = [b.wave_vac, r.wave_vac]
     flux = [b.cal_spec, r.cal_spec]
     ivar = [b.cal_ivar, r.cal_ivar]
  endelse 
  if arg_present(sensfunc) then sensfunc = [b.sensfunc, r.sensfunc]
  if arg_present(senswave) then senswave = [b.senswave, r.senswave]
  if arg_present(wave_air) then wave_air = [b.lambda, r.lambda]
end
