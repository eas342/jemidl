Function filter, name
  common jem$_filter, _filter_struct
  resolve_routine, 'spectrum__define', /no_recompile
  resolve_routine, 'dopplershift__define', /no_recompile
  if n_elements(_filter_struct) eq 0 then begin
     filterdir=getenv('JEM_REF')+'/FILTERS/'
     restore, filterdir+'filters.sav'
     resolve_obj, filters.landoltu
     _filter_struct = filters
  endif
  w = where( strupcase(name) eq tag_names(_filter_struct) )
  if w[0] eq -1 then return, 0
  return, { wave:_filter_struct.(w)->wavelength(), $
            throughput:_filter_struct.(w)->throughput() }
end
