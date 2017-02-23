Function BC03server::Init, $
   wavelengths=wavelengths, $
   sigmavbyc=sigmavbyc
  
  bc03_refdir = '/home/scpdata01/BC03/'
  file = file_search( djs_filepath( 'BC03.sav', root_dir=bc03_refdir ), count=ct )
  if (ct eq 0) then $
     message, 'Unable to find BC03.sav file'
  file = file[(reverse(sort(file)))[0]]
  restore, file

  self.age = ptr_new(_bc03spec_struct.age)
  self.metallicity = ptr_new(_bc03spec_struct.metallicity)
  

  if n_elements(wavelengths) eq 0 then $
     wavelengths = _bc03spec_struct.wave
  self.wave = ptr_new(wavelengths)

  age = _bc03spec_struct.age
  metallicity = _bc03spec_struct.metallicity

  sflux = fltarr( n_elements(wavelengths), $
                  n_elements(age), $
                  n_elements(metallicity) )
  
  tabinv, _bc03spec_struct.wave, wavelengths, waveindex
  if n_elements(sigmavbyc) eq 0 then begin
     for iage=0, n_elements(age)-1 do begin
        for imetal=0, n_elements(metallicity)-1 do begin
           sflux[*,iage,imetal] = interpolate( _bc03spec_struct.flux[*,iage,imetal], waveindex )
        endfor
     endfor
  endif else begin ;;this can't be right...  need dispersion in here somewhere...
     sigma = median(wavelengths)*sigmavbyc
     self.sigma = sigma
     kernel = findgen(11)-5
     kernel = exp(-kernel^2/(2*sigma^2))
     kernel /= total(kernel)
     for iage=0, n_elements(age)-1 do begin
        for imetal=0, n_elements(metallicity)-1 do begin
           flux = interpolate( _bc03spec_struct.flux[*,iage,imetal], waveindex )
           sflux[*,iage,imetal] = convol(flux, kernel, /edge_truncate)
        endfor
     endfor
  endelse
  self.flux = ptr_new(sflux)
  return, 1
end

Function BC03Server::Flux, metal=metal, age=age
  if n_elements(metal) eq 0 and n_elements(age) eq 0 then return, *self.flux
  if n_elements(metal) eq 0 and n_elements(age) ne 0 then begin
     flux = fltarr(n_elements(*self.wave), n_elements(*self.metallicity))
     tabinv, *self.age, age, ageindex
     for imetal=0, n_elements(*self.metallicity)-1 do begin
        flux[*,imetal] = interpolate((*self.flux)[*,*,imetal], ageindex)
     endfor
     return, flux
  endif
  if n_elements(metal) ne 0 and n_elements(age) eq 0 then begin
     flux = fltarr(n_elements(*self.wave), n_elements(*self.age))
     tabinv, *self.metallicity, metal, metalindex
     for iage=0, n_elements(*self.age)-1 do begin
        flux[*,iage] = interpolate((*self.flux)[*,iage,*], metalindex)
     endfor
     return, flux
  endif
  if n_elements(metal) ne 0 and n_elements(age) ne 0 then begin
     tabinv, *self.age, age, ageindex
     tabinv, *self.metallicity, metal, metalindex
     flux = interpolate(*self.flux, ageindex, metalindex)
     return, flux
  endif
end

Pro BC03Server::CleanUp
  ptr_free, self.flux, self.age, self.metallicity, self.wave
end

Pro BC03Server__Define, struct
  struct = { BC03Server, $
             Inherits JEMobject, $
             flux: ptr_new(), $
             age: ptr_new(), $
             metallicity: ptr_new(), $
             wave:ptr_new(), $
             sigma:0.0 }
end
