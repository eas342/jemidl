Function imag, spec, frame=frame
  common jem_$izmags, _izmag_struct
  resolve_routine, 'spectrum__define', /no_recompile
  resolve_routine, 'dopplershift__define', /no_recompile
  if n_elements(_izmag_struct) eq 0 then begin
     restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/data/filters.sav'
     restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/data/zeropointspecs.sav'
     resolve_obj, filters.acs_f775w
     resolve_obj, ABspec
     resolve_routine, 'dopplershift__define', /no_recompile
     izp = obj_new('zeropoint', ABspec, filters.acs_f775w)
     zzp = obj_new('zeropoint', ABspec, filters.acs_f850lp)
     imag = obj_new('magnitude', ABspec, filters.acs_f775w, izp)
     zmag = obj_new('magnitude', ABspec, filters.acs_f850lp, zzp)
     _izmag_struct = { ABspec:ABspec, $
                       ACS_F775W: filters.acs_f775w, $
                       ACS_F850LP: filters.acs_f850lp, $
                       IZP: izp, $
                       zzp: zzp, $
                       imag: imag, $
                       zmag: zmag }
  endif

  if n_elements(frame) eq 0 then frame='rest'

  junk = _izmag_struct.imag->spectrum( spectrum=spec )
  imag = _izmag_struct.imag->magnitude( frame=frame )
  return, imag
end
