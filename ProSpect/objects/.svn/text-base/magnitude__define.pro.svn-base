Function Magnitude::Init, $
   spectrum, $
   filter, $
   zeropoint

  self.zeropoint = zeropoint
  if not self->SynPhot::Init( spectrum, filter ) then return, 0
  return, 1
end

Function Magnitude::Magnitude, $
   _extra=extra

  zp_mag = self.zeropoint->base_mag( _extra=extra )
  spec_mag = self->base_mag( _extra=extra )
  return, spec_mag - zp_mag
end

Function Magnitude::Brightness_ratio, $
   brightness=brightness, $
   _extra=extra

  if n_elements( brightness ) eq 0 then $
     brightness = self->brightness( _extra=extra )
  zp_brightness = self.zeropoint->brightness( _extra=extra )
  return, brightness/zp_brightness
end

Function Magnitude::Brightness_ratio_ivar, $
   brightness=brightness, $
   brightness_ivar=brightness_ivar, $
   _extra=extra

  if n_elements(brightness) eq 0 then $
     brightness = self->brightness( _extra=extra )
  if n_elements(brightness_ivar) eq 0 then $
     brightness_ivar = self->brightness_ivar( _extra=extra )

  zp_brightness_ivar = self.zeropoint->brightness_ivar( _extra=extra )
  zp_brightness = self.zeropoint->brightness( _extra=extra )
  term1 = 1./(brightness_ivar*zp_brightness^2)
  term2 = 1./(zp_brighness_ivar)*brightess^(2.d)*zp_brightness^(-4.d)
  return, 1./(term1+term2)
end

Function Magnitude::Magnitude_errplus, $
   _extra=extra
  brightness = self->brightness( _extra=extra )
  brightness_ivar = self->brightness_ivar( _extra=extra )

  brightness_ratio = self->brightness_ratio( brightness=brightness, _extra=extra )
  brightness_ratio_ivar = self->brightness_ratio_ivar( brightness=brightness, $
                                                       brightness_ivar=brightness_ivar, $
                                                       _extra=extra )
  dbrightness_ratio = (brightness_ratio_ivar)^(-0.5)
  return, -2.5*alog10( 1. - dbrightness_ratio/brightness_ratio )
end

Function Magnitude::Magnitude_errminus, $
   _extra=extra
  brightness_ratio = self->brightness_ratio( _extra=extra )
  dbrightness_ratio = (self->brightness_ratio_ivar( _extra=extra ))^(-0.5)
  return, -2.5*alog10( 1. + dbrightness_ratio/brightness_ratio )
end

Pro Magnitude__Define
  struct = { Magnitude, $
             Inherits SynPhot, $
             zeropoint: Obj_New() }
end
