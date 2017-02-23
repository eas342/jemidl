Function fitgal, $
   i775, $
   z850, $
   zgal, $
   zcluster=zcluster, $
   ebv=ebv, $
   zform=zform, $
   age=age
  
  if n_elements(zcluster) eq 0 then zcluster=zgal
  if n_elements(ebv) eq 0 then ebv=0.d
  if n_elements(age) eq 0 then begin
     if n_elements(zform) eq 0 then age=2.5d9 $
     else age = galage( zcluster, zform, /silent )
  endif  

  metal = [0.0001d, 0.0004d, 0.004d, 0.008d, 0.02d, 0.05d]
  Eiz = dblarr( n_elements(metal) )
  
  lumratio = (lumdist( 1.0, /silent )/lumdist( zcluster, /silent ) )^2

  restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/filters.sav'
  restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/zeropointspecs.sav'
  resolve_routine, 'spectrum__define', /com, /no
  resolve_routine, 'filter__define', /com, /no
  resolve_routine, 'dopplershift__define', /com, /no

  i775AB_zp = obj_new( 'zeropoint', ABspec, filters.ACS_F775W )
  z850AB_zp = obj_new( 'zeropoint', ABspec, filters.ACS_F850LP )

  imag = obj_new( 'magnitude', ABspec, filters.ACS_F775W, i775AB_zp )
  zmag = obj_new( 'magnitude', ABspec, filters.ACS_F850LP, z850AB_zp )
  for i=0, n_elements(metal)-1 do begin
     spec = bc03spec( age, metal[i], z=zgal )
     ;;dust stuff here
     junk = imag->spectrum( spectrum=spec )
     junk = zmag->spectrum( spectrum=spec )
     Eiz[i] = (imag->base_mag(frame='obs')+25.678 - zmag->base_mag(frame='obs')-24.867) $
              - (i775-z850)
     obj_destroy, spec
  endfor
;  print, Eiz

  metal_out = interpol( metal, Eiz, 0., /spline )
  if metal_out lt min(metal) or metal_out gt max(metal) then begin
     return, {metal:-1.d, lum:-1.d}
  endif ;;panic
  
  spec = bc03spec( age, metal_out, z=zgal )
  ;; dust again
  junk = zmag->spectrum( spectrum=spec )
  zmag1 = zmag->magnitude(frame='obs')
  zmag1 -= 2.5*alog10(lumratio)
  lum = 10.^(0.4*(zmag1-z850))
  ;; return shit.
  obj_destroy, spec
  obj_destroy, imag, zmag
  obj_destroy, i775AB_zp, z850AB_zp
  return, {metal:metal_out, lum:lum}
end
