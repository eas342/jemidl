Pro scalespectra, spectra, _extra=extra
  for ispec=0, n_elements(spectra)-1 do begin
     tempspec = scalespec(spectra[ispec], refspec=spectra[0], _extra=extra)
     obj_destroy, spectra[ispec]
     spectra[ispec] = tempspec
  endfor
end
