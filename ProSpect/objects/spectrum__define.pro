Function Spectrum::Init, $
   wavelength, $
   flux, $
   ivar=ivar, $
   unit=unit, $
   orig_z=orig_z, $
   new_z=new_z

  if n_elements( unit ) eq 0 then unit='flambda'
  if n_elements( orig_z ) eq 0 then orig_z = Obj_New('DopplerShift', z=0.)
  if n_elements( new_z ) eq 0 then new_z = Obj_New('DopplerShift', z=0.)
  resolve_obj, orig_z

  ;1. Do consistency checks
  if n_elements( wavelength ) ne n_elements( flux ) then $
     message, 'Wavelength and flux arrays must have the same length!'
  if not member( unit, ['flambda','clambda'] ) then $
     message, 'Flux unit '+unit+' not recognized.'
  ;2. Sort by wavelength
  sort = sort( wavelength )
  self.wavelength = Ptr_New(wavelength[sort])
  self.flux = Ptr_New(self->convert_flux( unit, 'flambda', flux[sort], *self.wavelength ))
  ;3. Add redshift information
  self.new_z = new_z
  self.orig_z = orig_z
  ;4. De-redshift if orig_z ne 0.
  if self.orig_z->z() ne 0. then begin
     factor = self.orig_z->factor()
     *self.wavelength /= factor
     *self.flux *= factor^3
  endif
  ;5. Now take care of optional variance
  if n_elements(ivar) ne 0 then begin
     if n_elements(wavelength) ne n_elements(ivar) then $
        message, 'Wavelength and ivar arrays must have the same length!'
     stdev = ivar^(-0.5)
     self.ivar = Ptr_New((self->convert_flux( unit, 'flambda', stdev[sort], *self.wavelength))^(-2.d))
     if self.orig_z->z() ne 0. then begin
        factor = self.orig_z->factor()
        *self.ivar /= factor^6.d
     endif
  endif
  return, 1
end

Function Spectrum::Convert_Flux, $
   old_unit, $
   new_unit, $
   flux, $
   wavelength
  
h_erg_s = 6.626068e-27
c_AA_s = 2.99792e18

  if old_unit eq new_unit then $
     new_flux = flux $
  else if old_unit eq 'clambda' and new_unit eq 'flambda' then $
     new_flux = h_erg_s*c_AA_s*flux/wavelength $
  else if old_unit eq 'flambda' and new_unit eq 'clambda' then $
     new_flux = 1.d/(h_erg_s*c_AA_s)*flux*wavelength

  return, new_flux
end

Function Spectrum::Wavelength, $
   frame=frame

  if n_elements( frame ) eq 0 then frame = 'rest'
  wdata = *self.wavelength
  if frame eq 'obs' then begin
     factor = self.new_z->factor()
     wdata *= factor
  endif
  return, wdata
end

Function Spectrum::Flux, $
   unit=unit, $
   frame=frame
  
  if n_elements( unit ) eq 0 then unit = 'flambda'
  if n_elements( frame ) eq 0 then frame = 'rest'

  fdata = *self.flux
  wdata = *self.wavelength
  if frame eq 'obs' then begin
     factor = self.new_z->factor()
     fdata /= factor^3
     wdata *= factor
  endif

  f_final = self->convert_flux( 'flambda', unit, fdata, wdata )
  return, f_final
end

Function Spectrum::Ivar, $
   unit=unit, $
   frame=frame

  if not ptr_valid(self.ivar) then $
     message, 'Inverse variance array not defined!'

  if n_elements( unit ) eq 0 then unit = 'flambda'
  if n_elements( frame ) eq 0 then frame = 'rest'

  idata = *self.ivar
  wdata = *self.wavelength
  sdata = idata^(-0.5)
  if frame eq 'obs' then begin
     factor = self.new_z->factor()
     sdata /= factor^3
     wdata *= factor
  endif
  s_final = self->convert_flux( 'flambda', unit, sdata, wdata )
  i_final = s_final^(-2.d)
  return, i_final
end     

Function Spectrum::Orig_z
  return, self.orig_z
end

Function Spectrum::New_z, $
   new_z=new_z

  if n_elements(new_z) ne 0 then self.new_z=Obj_New('DopplerShift',new_z)
  return, self.new_z
end

Pro Spectrum::CleanUp
  ptr_free, self.wavelength
  ptr_free, self.flux
  ptr_free, self.ivar
  obj_destroy, self.orig_z
  obj_destroy, self.new_z
end

Pro Spectrum__Define
  struct = { Spectrum, $
             wavelength: Ptr_New(), $
             flux: Ptr_New(), $
             ivar: Ptr_New(), $
             orig_z: Obj_New(), $
             new_z: Obj_New() }
end
