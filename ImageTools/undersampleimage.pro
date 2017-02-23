;+
; NAME:
;	UNDERSAMPLEIMAGE
;
; PURPOSE:
;       This function returns the undersampling of an image.
; CALLING SEQUENCE:
;	Result = UNDERSAMPLEIMAGE(image, undersample_factor)
; INPUTS:
;       image: 2d array to be undersampled
;       undersample_factor: factor by which to undersample the image
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       Result: The undersampled image.
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       JEM, Jan, 2009. Written
;-

Function undersampleimage, image, undersample_factor
  sx = ((size(image,/dim))[0])/undersample_factor
  sy = ((size(image,/dim))[1])/undersample_factor
  undersampleimage = dblarr(sx, sy)
  for ix=0, sx-1 do begin
     for iy=0, sy-1 do begin
        undersampleimage[ix,iy]=mean(image[undersample_factor*(ix):undersample_factor*(ix+1)-1, $
                                     undersample_factor*(iy):undersample_factor*(iy+1)-1], /nan)
     endfor
  endfor
  return, undersampleimage
end
