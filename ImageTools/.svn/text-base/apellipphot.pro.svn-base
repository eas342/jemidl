;+
; NAME:
;	APELLIPPHOT
;
; PURPOSE:
;       Perform elliptical aperture photometry on an image.
; CALLING SEQUENCE:
;	Result = APELLIPPHOT(image, ...)
; INPUTS:
;       image: An array of pixels...
;       x/y: the center of the ellipse in josh coordinates
;       boa: the minor axis to major axis ratio or the ellipse
;       pa: the position angle of the ellipse
;       radius: vector of major axis values within which to measure
;               the flux
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       Result: The oversampled sum of all pixels within the specified
;       ellipses...  (not in magnitudes!, just flux!)
; EXAMPLE:
;
; NOTES:
;       Josh coords place the origin at the lower left corner of the
;       lower left pixel.
;
; MODIFICATION HISTORY:
;       JEM, Jan, 2009. Written
;-

Function ApEllipPhot, image, x, y, boa, pa, radius
  oversample_factor=5
  overimage = oversampleimage(image, oversample_factor)
  dist_ellipse, dist, size(overimage,/dim), $
                x*oversample_factor-0.5, y*oversample_factor-0.5, $
                boa, pa-90., /double
  flux = dblarr(n_elements(radius))*!values.d_nan
  for ir=0, n_elements(radius)-1 do begin
     if ~finite(radius[ir]) then continue
     mask = dist le radius[ir]*oversample_factor
     flux[ir] = total(overimage*mask, /nan)
     flux[ir] *= !dpi*radius[ir]*radius[ir]*boa/total(mask)
;     w=where(dist le radius[ir]*oversample_factor)
;     if w[0] eq -1 then flux[ir] = !values.f_nan $
;     else begin
;        flux[ir] = total(overimage[w], /nan)
;        flux[ir] *= 3.14159*radius[ir]*radius[ir]*boa/n_elements(w)
;     endelse
  endfor
  return, flux
end
