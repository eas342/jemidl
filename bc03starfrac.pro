Function BC03StarFrac, $
   age, $
   metal

  ;; store _BC03spec_struct for future use after opening...
  common jem_$bc03starfrac, _bc03spec_struct
  if n_elements(_bc03spec_struct) eq 0 then begin
     bc03_refdir = '/home/scpdata01/BC03/'
     file = file_search(djs_filepath('BC03.sav', root_dir=bc03_refdir), count=ct)
     if (ct eq 0) then $
        message, 'Unable to find BC03.sav file'
     file = file[(reverse(sort(file)))[0]]
     restore, file
  endif
  
  ;; find the effective indices for the age and metallicities...
  tabinv, _bc03spec_struct.age, age, ageindex
  tabinv, _bc03spec_struct.metallicity, metal, metalindex
  ;; now do the bilinear interpolation...
  starfrac = interpolate( _bc03spec_struct.starmass, ageindex, metalindex )
  return, starfrac
end
