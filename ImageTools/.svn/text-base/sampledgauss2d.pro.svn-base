;+
; NAME:
;	SAMPLEDGAUSS2D
;
; PURPOSE:
;       This function makes an image of a 2d gaussian factoring in the
;       effects of sampling.
;       
; CALLING SEQUENCE:
;	Result = SAMPLEDGAUSS2d(xs, ys, param)
; INPUTS:
;       xs/ys: coords of pixel centers in resulting image.  each is a 1d vector.
;       param: an array with values as follows:
;              param[0] : offset
;              param[1] : amplitude
;              param[2] : x-sigma
;              param[3] : y-sigma
;              param[4] : x-center
;              param[5] : y-center
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       Result: The properly sampled image of a gaussian
; EXAMPLE:
;
; NOTES:
;       Josh coords place the origin at the lower left corner of the
;       lower left pixel.
;
; MODIFICATION HISTORY:
;       JEM, Jan, 2009. Written
;-

Function SampledGauss2d, xs, ys, param
  oversample_factor = 9
  halfbin = (oversample_factor-1)/2
  
  newx = rebin(xs, n_elements(xs),oversample_factor,/sample)
  offx = rebin(transpose((findgen(oversample_factor)-halfbin)/oversample_factor),n_elements(xs),oversample_factor,/sample)
  newxs = reform(transpose(newx+offx), n_elements(xs)*oversample_factor)

  newy = rebin(ys, n_elements(ys),oversample_factor,/sample)
  offy = rebin(transpose((findgen(oversample_factor)-halfbin)/oversample_factor),n_elements(ys),oversample_factor,/sample)
  newys = reform(transpose(newy+offy), n_elements(ys)*oversample_factor)

  make_2d, newxs, newys
  u = ((newxs-param[4])/param[2])^2+((newys-param[5])/param[3])^2
  overmodel = param[0]+param[1]*exp(-0.5*u)
  model = undersampleimage(overmodel,oversample_factor)
  return, model
end
