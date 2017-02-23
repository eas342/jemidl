;+
; NAME:
;	MEMBER
;
; PURPOSE:
;       This function determines whether each of the items in ELTS is
;       a member of SET or not.
; CALLING SEQUENCE:
;	Result = MEMBER(ELTS, SET)
; INPUTS:
;       elts: A arrat of test elements.
;       set: A set.
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       Result: Boolean array with the results of the member test.
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       JEM, July, 2008. Written
;-

Function Member, elts, set
  result = bytarr(n_elements(elts))
  for j=0L, n_elements(elts)-1 do result[j] = (total(elts[j] eq set) ge 1)
  return, result
end
