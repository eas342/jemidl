;+
; NAME:
;	RM_ELT
;
; PURPOSE:
;       This function removes specified elements from an array.
; CALLING SEQUENCE:
;	Result = RM_ELT(array, wbad)
; INPUTS:
;       array: Array from which elements are to be removed.
;       wbad: Indices of elements to be removed
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       Result: Array with specified elements removed
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       JEM, May, 2010. Written
;-

Function Rm_elt, array, wbad
  wgood=indgen(n_elements(array))
  junk=where(member(wgood, wbad), complement=w)
  return, array[w]
end
