Function ModeWait, a, j

;  Function to estimate the mode of an array.
;	Uses a method described in Press, et al, Numerical Recipes, p. 463.
;
; INPUT:
;	a = array for which we want estimate of mode.
;	j = window size to slide across sorted data.  
;		(This is difficult to explain without reading the section
;		of Numerical Recipes.  Then it should be obvious.)
; OUTPUT:
;	returns estimate of mode.
;


na = n_elements(a)
if n_elements(j) ne 1 then j = na/10
j = j>3

sorted = a[sort(a)]
usemin = min( (sorted[j:*] - sorted[0:na-1-j]), minhere)
modevalue = (sorted[minhere] + sorted[minhere+j])/2.

return, modevalue
end
