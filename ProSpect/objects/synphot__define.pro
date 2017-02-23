;+
; NAME:
;	CLASS SYNPHOT
;
; PURPOSE:
;       Implements SynPhot class.
;
; DATA:
;       spectrum : Spectrum on which to perform synthetic photometry.
;       filter   : Filter with which to perform synthetic photometry.
;       
; METHODS:
;       Init(spectrum, filter) : Initialize with given spectrum and filter.
;       SynPhot : Return product of filter throughput and spectrum as a
;                 function of wavelength.
;       Ivar    : Return inverse variance of SynPhot as a function of wavelength.
;       Brightness : Return integral of SynPhot
;       Brightness_Ivar : return inverse variance of brightness
;       Base_mag : Return -2.5 * log10( brightness )
;       Base_magerrplus/minus : Return asymmetric error bars of base_mag
;       
;       Accesors:
;             Spectrum : Set or return spectrum.
;             Filter   : Set or return filter.
;
; MODIFICATION HISTORY:
;       JEM, Aug, 2008. Written
;-


Function SynPhot::Init, $
   spectrum, $
   filter

  resolve_obj, spectrum
  resolve_obj, filter
  self.spectrum = spectrum
  self.filter = filter

  return, 1
end

Function SynPhot::SynPhot, $
   frame=frame, $
   unit=unit
  
  if n_elements( frame ) eq 0 then frame = 'obs'
  if n_elements( unit ) eq 0 then unit = 'clambda' ;e.g, clambda for CCDs, flambda for phototubes...

  sp_flux = self.spectrum->flux( unit=unit, frame=frame )
  sp_wave = self.spectrum->wavelength( frame=frame )
  sp_throughput = self.filter->throughput( wavelength=sp_wave )
  sp_values = sp_flux*sp_throughput
  return, sp_values
end

Function SynPhot::Ivar, $
   frame=frame, $
   unit=unit

  if n_elements( frame ) eq 0 then frame = 'obs'
  if n_elements( unit ) eq 0 then unit = 'clambda'

  sp_ivar = self.spectrum->ivar( unit=unit, frame=frame )
  sp_wave = self.spectrum->wavelength( frame=frame )
  sp_throughput = self.filter->throughput( wavelength=sp_wave )
  sp_ivar_values = sp_ivar*sp_throughput^(-2.d)
  return, sp_ivar_values
end

Function SynPhot::Brightness, $
   frame=frame, $
   unit=unit
  
  if n_elements( frame ) eq 0 then frame = 'obs'
  if n_elements( unit ) eq 0 then unit = 'clambda'

  ;check that filter fits within spectral range
  sw = self.spectrum->wavelength( frame=frame )
  fw = self.filter->wavelength()
  if min( fw ) lt min( sw ) or max( fw ) gt max( sw ) then $
     message, 'Filter not within bounds of Spectrum wavelength coverage.'
  sp_wave = sw
  sp_dwave = (sp_wave - shift(sp_wave, 1))[1:*]
  sp = (self->synphot( frame=frame, unit=unit ))[1:*]
  brightness = total(sp*sp_dwave)
  return, brightness
end

Function SynPhot::Brightness_ivar, $
   frame=frame, $
   unit=unit

  if n_elements( frame ) eq 0 then frame = 'obs'
  if n_elements( unit ) eq 0 then unit = 'clambda'

  ;check that filter fits within spectral range
  sw = self.spectrum->wavelength( frame=frame )
  fw = self.filter->wavelength()
  if min( fw ) lt min( sw ) or max( fw ) gt max( sw ) then $
     message, 'Filter not within bounds of Spectrum wavelength coverage.'
  sp_dwave = (sw - shift(sw, 1))[1:*]
  sp_ivar = (self->ivar( frame=frame, unit=unit ))[1:*]
  brightness_ivar = 1./(total( (1.d/sp_ivar)*sp_dwave^(2.d) ))
  return, brightness_ivar
end

Function SynPhot::Base_mag, $
   frame=frame, $
   unit=unit

  if n_elements( frame ) eq 0 then frame = 'rest'
  if n_elements( unit ) eq 0 then unit = 'clambda'
  return, -2.5*alog10( self->brightness( frame=frame, unit=unit ) )
end

Function SynPhot::Base_magerrplus, $
   frame=frame, $
   unit=unit

  if n_elements( frame ) eq 0 then frame = 'rest'
  if n_elements( unit ) eq 0 then unit = 'clambda'
  
  brightness = self->brightness( frame=frame, unit=unit )
  dbrightness = (self->brightness_ivar( frame=frame, unit=unit ))^(-0.5d)
  
  return, -2.5*alog10( 1. - dbrightness/brightness )
end

Function SynPhot::Base_magerrminus, $
   frame=frame, $
   unit=unit

  if n_elements( frame ) eq 0 then frame = 'rest'
  if n_elements( unit ) eq 0 then unit = 'clambda'
  
  brightness = self->brightness( frame=frame, unit=unit )
  dbrightness = (self->brightness_ivar( frame=frame, unit=unit ))^(-0.5d)
  
  return, -2.5*alog10( 1. + dbrightness/brightness )
end

Function SynPhot::Spectrum, $
   spectrum=spectrum

  if n_elements(spectrum) ne 0 then self.spectrum=spectrum
  return, self.spectrum
end

Function SynPhot::Filter, $
   filter=filter

  if n_elements(filter) ne 0 then self.filter=filter
  return, self.filter
end

Pro SynPhot__Define
  struct = { SynPhot, $
             spectrum: Obj_New(), $
             filter: Obj_New() }
end
