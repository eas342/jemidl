;+
; NAME:
;	OVERSAMPLEIMAGE
;
; PURPOSE:
;       This function oversamples an image.
; CALLING SEQUENCE:
;	Result = OVERSAMPLEIMAGE(image, oversample_factor, [binfrac=binfrac, interp=interp])
; INPUTS:
;       image: 2d array to be oversampled
;       oversample_factor: factor by which to oversample the image
; KEYWORD PARAMETERS:
;       binfrac:  Set this keyword to set the edge subpixels of
;                 each original pixel to NaN
;       interp:   Set this to use cubic interpolation feature of the
;                 IDL function congrid.
; OUTPUTS:
;       Result: The oversampled image
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       JEM, Jan, 2009. Written
;-
Function OverSampleImage, image, oversample_factor, binfrac=binfrac, interp=interp
  sz=size( image, /dim )
  xsize = sz[0]
  ysize = sz[1]
  if keyword_set(interp) then out = congrid( image, oversample_factor*xsize, $
                                             oversample_factor*ysize, /center, cubic=-0.5 ) $
  else $
     out = rebin( image, oversample_factor*xsize, oversample_factor*ysize, /sample )

  if keyword_set(binfrac) then begin
     for i=0,xsize-1 do begin
        column1 = i*oversample_factor
        column2 = (i+1)*oversample_factor-1
        out[[column1,column2],*]=!values.f_nan
     endfor
     for i=0,ysize-1 do begin
        row1 = i*oversample_factor
        row2 = (i+1)*oversample_factor-1
        out[*,[row1,row2],*]=!values.f_nan
     endfor
  endif
  return, out
end
