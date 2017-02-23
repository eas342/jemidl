;+
; NAME:
;	CLASS DOPPLERSHIFT
;
; PURPOSE:
;       Implements DopplerShift class.  Use for easy conversions between
;       multiplier, redshift, and velocity.
;
; DATA:
;       Z : store the redshift.
; METHODS:
;       Init(factor, z, v) : Initialize with one of factor, redshift, or velocity.
;       Accesors: 
;             Factor : Set or return factor
;             Z      : Set or return redshift
;             V      : Set or return velocity
;
; MODIFICATION HISTORY:
;       JEM, Aug, 2008. Written
;-


Function DopplerShift::Init, $
   factor=factor, $
   z=z, $
   v=v

  c_kms = 299792.d
  if n_elements(factor) ne 0 then self.z = factor - 1.d
  if n_elements(v) ne 0 then self.z = sqrt((1.d + v/c_kms)/(1.d - v/c_kms))
  if n_elements(z) ne 0 then self.z = z
  return, 1
end

Function DopplerShift::Factor, $
   factor=factor

  if n_elements(factor) eq 1 then self.z = factor - 1.d
  return, self.z + 1
end

Function DopplerShift::Z, $
   z=z

  if n_elements(z) eq 1 then self.z = z
  return, self.z
end

Function DopplerShift::V, $
   v=v

  c_kms = 299792.d
  if n_elements(v) eq 1 then self.z = sqrt((1.d + v/c_kms)/(1.d - v/c_kms))
  return, c_kms * ((1.d + self.z)^2 - 1.d)/((1.d + self.z)^2 + 1.d)
end

Pro DopplerShift__Define
  struct = { DopplerShift, $
             z: 0.d }
end
