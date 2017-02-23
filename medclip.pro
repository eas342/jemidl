;; Routine to identify and replace large pixel outliers (wrt local variance)

Function medclip, spec, window, nsig
  medspec = medsmooth(spec, window) ;; simple median filtered spectrum
  madspec = medsmooth(abs(spec-medspec), 2*window) ;; estimate of local variance
  madspec = madspec > median(madspec) ;; need to protect against zero local variance
  out = spec
  w=where(abs((spec-medspec)/madspec) gt nsig) ;; replace large outliers with the medspec
  if w[0] ne -1 then out[w] = medspec[w]
  return, out
end
