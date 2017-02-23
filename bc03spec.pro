;+
; NAME:
;	BC03SPEC
;
; PURPOSE:
;       This function returns the SED of a simple stellar population from the
;       Bruzual/Charlot models.  The return type is a ProSpect spectrum object.
; CALLING SEQUENCE:
;	spec = BC03Spec( age, metal, z=z )
; INPUTS:
;       age: The age in yr of the desired model.
;       metal: The metallicity Z, by fractional mass (solar = 0.0177), of the desired model.  
; KEYWORD PARAMETERS:
;       z: Optionally redshift the spectrum.
;       erg: Optionally convert spectrum to ergs instead of L_sun
;       origmass: Optionally scale spectrum to galaxy of original
;       stellar mass origmass
; OUTPUTS:
;       Result: a ProSpect "spectrum" class object.
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       JEM, July, 2008. Written
;-

Function BC03Spec, $
   age, $
   metal, $
   z=z, $
   erg=erg, $
   origmass=origmass

  resolve_routine, 'spectrum__define', /no_recompile
  resolve_routine, 'dopplershift__define', /no_recompile

  ;; store _BC03spec_struct for future use after opening...
  common jem_$bc03spec, _bc03spec_struct
  if n_elements(z) eq 0 then z=0.d
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
  flux=interpolate( _bc03spec_struct.flux, ageindex, metalindex )
  if n_elements(erg) ne 0 then flux *= 3.826d33
  if n_elements(origmass) ne 0 then flux *= origmass
  return, obj_new( 'spectrum', _bc03spec_struct.wave, flux, $
                   new_z=obj_new( 'DopplerShift', z=z ), unit='flambda' )
end
