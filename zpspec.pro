Function zpspec, name
  common jem$_zpspec, _zpspec_struct
  ;; resolve_routine, 'spectrum__define', /no_recompile
  ;; resolve_routine, 'dopplershift__define', /no_recompile
  if n_elements(_zpspec_struct) eq 0 then begin
     zpdir=getenv('JEM_REF')+'/SPECTRA/'
     readcol, zpdir+'vega.dat', vwave, vflux, /silent
     readcol, zpdir+'AB.dat', abwave, abflux, /silent
     ;; restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/data/zeropointspecs.sav'
     ;; resolve_obj, ABspec
     ;; resolve_routine, 'dopplershift__define', /no_recompile
     ;; _zpspec_struct = { ABspec:ABspec, $
     ;;                    vegaspec:vegaspec, $
     ;;                    STspec:STspec }
     _zpspec_struct = { AB:ptr_new({wave:abwave, flux:abflux}), $
                        vega:ptr_new({wave:vwave, flux:vflux}) }
  endif
  w = where( strupcase(name) eq tag_names(_zpspec_struct) )
  if w[0] eq -1 then return, 0
  ;; return, { wave:_zpspec_struct.(w)->wavelength(), $
  ;;           flux:_zpspec_struct.(w)->flux() }
  return, { wave:(*_zpspec_struct.(w)).wave, $
            flux:(*_zpspec_struct.(w)).flux }

end
