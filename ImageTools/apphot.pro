;+
; NAME:
;	APPHOT
;
; PURPOSE:
;       Perform circular aperture photometry on an image.
; CALLING SEQUENCE:
;	Result = APPHOT(image, ...)
; INPUTS:
;       image: An array of pixels...
;       x/y: the center of the ellipse in josh coordinates
;       radius: vector of radii within which to measure
;               the flux
;      
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       Result: The oversampled sum of all pixels within the specified
;       circles...  (not in magnitudes!, just flux!)
; EXAMPLE:
;
; NOTES:
;       Josh coords place the origin at the lower left corner of the
;       lower left pixel.
;
; MODIFICATION HISTORY:
;       JEM, Jan, 2009. Written
;-

Function ApPhot, image, x, y, radius
;  oversample_factor=5
;  overimage = oversampleimage(image, oversample_factor)
;  dist_ellipse, dist, size(overimage,/dim), $
;                x*oversample_factor-0.5, y*oversample_factor-0.5, $
;                1, 0, /double
;  flux = dblarr(n_elements(radius))
;  for ir=0, n_elements(radius)-1 do begin
;     w=where(dist LE radius[ir]*oversample_factor)
;     if w[0] eq -1 then flux[ir] = !values.f_nan $
;     else begin 
;        flux[ir] = total(overimage[w])
;        flux[ir] *= !dpi*radius[ir]*radius[ir]/n_elements(w)
;     endelse
;  endfor
;  return, flux
  return, ApEllipPhot( image, x, y, 1., 0., radius )
end
