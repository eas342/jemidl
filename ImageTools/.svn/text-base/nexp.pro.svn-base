;+
; NAME:
;	NEXP
;
; PURPOSE:
;       Convert a multidrizzle context image to a number of exposures image.
; CALLING SEQUENCE:
;	Result = NEXP(data, ...)
; INPUTS:
;       context: A multidrizzle context image.
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       Result: The number of exposures in each pixel.
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       JEM, Jan, 2009. Written
;-


Function NExp, context
  if size( context, /n_dim ) eq 2 then nslice=1 $
  else nslice = (size( context, /dim ) )[2]
  nx = (size( context, /dim ))[0]
  ny = (size( context, /dim ))[1]
  nexp = intarr( nx, ny )
  for islice=0, nslice-1 do begin
     for i=0, 31 do begin
        nexp += (context[*,*,islice] and 2L^i) eq 2L^i
     endfor
  endfor
  return, nexp
end
