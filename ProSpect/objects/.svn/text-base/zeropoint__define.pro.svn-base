Function ZeroPoint::Init, $
   spectrum, $
   filter, $
   unit=unit

  if n_elements(unit) eq 0 then unit='clambda'

  if ~self->SynPhot::Init( spectrum, filter ) then return, 0
  self.unit = unit
  return, 1
end

Function ZeroPoint::Base_mag, $
   force=force, $
   unit=unit, $
   _extra=extra ;needed because SynPhot::Base_mag has more keywords...

  if n_elements(unit) EQ 0 then unit = self.unit $
  else if unit ne self.unit then begin
     self->reset
     self.unit = unit
  endif
  if keyword_set(force) or ~self.mag_done then begin
     self.stored_mag = self->SynPhot::base_mag( unit=unit )
     self.mag_done = 1
  endif
  return, self.stored_mag
end

Function ZeroPoint::Base_magerrplus, $
   force=force, $
   unit=unit, $
   _extra=extra

  if n_elements(unit) EQ 0 then unit = self.unit $
  else if unit ne self.unit then begin
     self->reset
     self.unit = unit
  endif
  if keyword_set(force) or ~self.magerrplus_done then begin
     self.stored_magerrplus = self->SynPhot::base_magerrplus( unit=unit )
     self.magerrplus_done = 1
  endif
  return, self.stored_magerrplus
end

Function ZeroPoint::Base_magerrminus, $
   force=force, $
   unit=unit, $
   _extra=extra
   

  if n_elements(unit) EQ 0 then unit = self.unit $
  else if unit ne self.unit then begin
     self->reset
     self.unit = unit
  endif
  if keyword_set(force) or ~self.magerrminus_done then begin
     self.stored_magerrminus = self->SynPhot::base_magerrminus( unit=unit )
     self.magerrminus_done = 1
  endif
  return, self.stored_magerrminus
end

Function ZeroPoint::Brightness, $
   force=force, $
   unit=unit, $
   _extra=extra

  if n_elements(unit) eq 0 then unit = self.unit $
  else if unit ne self.unit then begin
     self->reset
     self.unit = unit
  endif
  if keyword_set(force) or ~self.brightness_done then begin
     self.stored_brightness = self->SynPhot::brightness( unit=unit )
     self.brightness_done = 1
  endif
  return, self.stored_brightness
end

Function ZeroPoint::Brightness_ivar, $
   force=force, $
   unit=unit, $
   _extra=extra
   

  if n_elements(unit) eq 0 then unit = self.unit $
  else if unit ne self.unit then begin
     self->reset
     self.unit = unit
  endif
  if keyword_set(force) or ~self.brightness_ivar_done then begin
     self.stored_brightness_ivar = self->SynPhot::brightness_ivar( unit=unit )
     self.brightness_ivar_done = 1
  endif
  return, self.stored_brightness_ivar
end


Pro ZeroPoint::Reset
  self.mag_done = 0
  self.magerrplus_done = 0
  self.magerrminus_done = 0
  self.brightness_done = 0
  self.brightnessivar_done = 0
end

Pro ZeroPoint__Define
  struct = { ZeroPoint, $
             Inherits SynPhot, $
             stored_mag:0.d, $
             stored_magerrplus:0.d, $
             stored_magerrminus:0.d, $
             stored_brightness:0.d, $
             stored_brightness_ivar:0.d, $
             mag_done: 0, $
             magerrplus_done: 0, $
             magerrminus_done: 0, $
             brightness_done: 0, $
             brightness_ivar_done: 0, $
             unit: '' }
end
