Function flag_o2_telluric_pixels, $
   geoWave ;; geocentric frame wavelengths

  flag = bytarr(n_elements(geoWave))
  flag or= geoWave gt 6270.2 and geoWave lt 6331.7
  flag or= geoWave gt 6862.1 and geoWave lt 6964.6 ;; B-band
;  flag or= geoWave gt 7143.3 and geoWave lt 7398.2
  flag or= geoWave gt 7585.8 and geoWave lt 7703.0 ;; A-band
;  flag or= geoWave gt 7887.6 and geoWave lt 8045.8
;  flag or= geoWave gt 8916.0 and geoWave lt 9929.8

  return, flag
end
