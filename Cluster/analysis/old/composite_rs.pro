Function composite_RS, specgals, metallum, zcluster, vel_disp=vel_disp, zform=zform, age=age
  if n_elements(age) eq 0 then begin
     if n_elements(zform) eq 0 then age=2.5d9 $
     else age = galage( zcluster, zform, /silent )
  endif  
  if n_elements(vel_disp) eq 0 then vel_disp=0. ;;km/s
  c_kms=299792.d

  restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/filters.sav'
  restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/zeropointspecs.sav'
  resolve_routine, 'spectrum__define', /com, /no
  resolve_routine, 'filter__define', /com, /no
  resolve_routine, 'dopplershift__define', /com, /no

  i775AB_zp = obj_new( 'zeropoint', ABspec, filters.ACS_F775W )
  z850AB_zp = obj_new( 'zeropoint', ABspec, filters.ACS_F850LP )
  imag = obj_new( 'magnitude', ABspec, filters.ACS_F775W, i775AB_zp )
  zmag = obj_new( 'magnitude', ABspec, filters.ACS_F850LP, z850AB_zp )

  lumratio = (lumdist( 1.0, /silent )/lumdist( zcluster, /silent ))^2

  for igal=0, n_elements(specgals)-1 do begin
     if metallum[igal].lum eq -1 then iz1 = {i775:0.d, z850:0.d} $
     else begin
        zgal = zcluster + specgals[igal].velocity/specgals[igal].vel_disp*vel_disp*(1.+zcluster)/c_kms
        spec = bc03spec( age, metallum[igal].metal, z=zgal )
        junk = imag->spectrum( spectrum = spec )
        junk = zmag->spectrum( spectrum = spec )
        i775 = imag->magnitude( frame='obs' )
        z850 = zmag->magnitude( frame='obs' )

        i775 -= 2.5*alog10(lumratio)
        z850 -= 2.5*alog10(lumratio)

        i775 -= 2.5*alog10(metallum[igal].lum)
        z850 -= 2.5*alog10(metallum[igal].lum)
        
        iz1 = {i775:i775, z850:z850}
     endelse
     if n_elements(iz) eq 0 then iz=iz1 $
     else iz=[iz,iz1]
  endfor
  return, iz
end
