;+
; NAME:
;	ERRIM
;
; PURPOSE:
;       Calculates from first principles an error image.
; CALLING SEQUENCE:
;	result = ERRIM(...)
; INPUTS:
;       image: An array of pixel values in [counts/sec]
;       time: An array of exposure times per each pixel. [sec]
;       gain: The average gain of all the pixels. [counts/e-]
;       rdnoise: The average read noise of each amplifier. [cts]
;       sky: The average sky value for each pixel. [cts/sec]
;       dark: The average dark current for each pixel.  [cts/sec]
;       nexp: An array of exposure counts for each pixel. 
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       Result: The error image. [cts/sec]
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       JEM, Jan, 2009. Written
;-

Function ErrIm, image, time, gain, rdnoise, sky, dark, nexp
  return, sqrt( image*time*gain + rdnoise*rdnoise*nexp + (sky+dark)*time*gain ) / (gain * time)
end
