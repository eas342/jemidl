Function flag_emission_line_pixels, $
   restWave ;; Rest-frame wavelength

  c = 299792.458d ;; speed of light in km/s

  lines = [4862.677d, 4960.295d, 5008.239d]              ;; Hbeta, [OIII]
  lines = [lines, 6564.614d, 4341.676d, 4102.884d]       ;; Ha, Hg, Hd
  lines = [lines, 3971.188d, 3890.143d]                  ;; H5, H6
  lines = [lines, 3727.092d, 3729.875d]                  ;; [OII]
  lines = [lines, 6549.86d, 6585.27d, 6718.29, 6732.68d] ;; [NII], [SII]

;; width/2 of masked gas emission region in km/s
  dv = [800d, 500d, 800d]            ;; Hbeta, [OIII]
  dv = [dv, 500d, 500d, 500d]        ;; Ha, Hg, Hd
  dv = [dv, 500d, 500d]              ;; H5, H6
  dv = [dv, 500d, 500d]              ;; [OII]
  dv = [dv, 500d, 500d, 500d, 500d]  ;; [NII], [SII]

  flag = bytarr(n_elements(restWave))

  for j=0, n_elements(lines)-1 do $
     flag or= restWave gt lines[j]*(1d - dv[j]/c) $
     and restWave lt lines[j]*(1d + dv[j]/c)

  return, flag
end
