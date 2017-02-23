Function coloragemetal, age, metal, z=z
  resolve_routine, 'spectrum__define', /com, /no
  resolve_routine, 'filter__define', /com, /no
  resolve_routine, 'dopplershift__define', /com, /no

  restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/data/filters.sav'
  restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/data/zeropointspecs.sav'

  if n_elements(z) eq 0 then z=0.
  izp = obj_new('zeropoint', ABspec, filters.acs_f775w)
  zzp = obj_new('zeropoint', ABspec, filters.acs_f850lp)
  spec = bc03spec(age, metal, z=z)
  imag = obj_new('magnitude', spec, filters.acs_f775w, izp )
  zmag = obj_new('magnitude', spec, filters.acs_f850lp, zzp )
  junk = imag->spectrum(spectrum=spec)
  junk = zmag->spectrum(spectrum=spec)
  color = imag->magnitude(frame='obs')-zmag->magnitude(frame='obs')
  obj_destroy, [spec, imag, zmag]
  obj_destroy, [izp, zzp, vegaspec, stspec, abspec]
  return, color
end
