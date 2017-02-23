;+
; NAME:
;	RR
;
; PURPOSE:
;       This function updates once any stopped widget programs.  Note that you 
;       need to type Ctrl-C and then execute 'return' to return to the execution
;       point of a stopped program when rr was started.
; CALLING SEQUENCE:
;	RR
; INPUTS:
; KEYWORD PARAMETERS:
; OUTPUTS:
; EXAMPLE:
; MODIFICATION HISTORY:
;       JEM, July, 2008. Written
;-
Pro Rr
  while 1 do begin & wait, 0.1 & void = widget_event(/NOWAIT) & endwhile
end
