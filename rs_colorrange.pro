Function RS_ColorRange, z, zform=zform, age=age1
  resolve_routine, 'spectrum__define', /com
  resolve_routine, 'filter__define', /com
  resolve_routine, 'dopplershift__define', /com

  restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/data/filters.sav'
  restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/data/zeropointspecs.sav'

  if n_elements(zform) eq 0 then zform=5.0
  age = galage(z,zform,/silent)
  if n_elements(age1) ne 0 then age=age1
  metals = [0.004, 0.008, 0.020, 0.050]
  izp = obj_new('zeropoint', vegaspec, filters.acs_f775w)
  zzp = obj_new('zeropoint', vegaspec, filters.acs_f850lp)
  color = fltarr(n_elements(metals))
  for imetal=0, n_elements(metals)-1 do begin
     spec = bc03spec(age, metals[imetal], z=z)
     imag = obj_new('magnitude', spec, filters.acs_f775w, izp )
     zmag = obj_new('magnitude', spec, filters.acs_f850lp, zzp )
     color[imetal] = imag->magnitude(frame='obs')-zmag->magnitude(frame='obs')
     obj_destroy, [spec, imag, zmag]
  endfor
  obj_destroy, [izp, zzp, vegaspec, stspec, abspec]
  return, minmax(color)
end
