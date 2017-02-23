; Generate Chebyshev polynomials

Function Cheb_coeffs, n
; Use recurrence relations to generate Chebyshev polynomial coefficients.
  if (n eq 0) then begin
    c = 1
    return, c
  endif
  c = intarr(n+1, n+1)
  c[0, 0] = 1
  c[1, 1] = 1
  for i = 2, n do begin
     c[i, *] = 2*shift(c[i-1,*], 1)-c[i-2, *]
  endfor
  return, c
end

;+
; NAME:
;	CHEB_ARRAY
;
; PURPOSE:
;       This function returns as an array the first n+1 Chebyshev
;       polynomials defined on the domain [-1,1] rescaled to the domain
;       of x.
; CALLING SEQUENCE:
;	Result = CHEB_ARRAY(n, x)
; INPUTS:
;       n: The degree of the highest Chebyshev polynomial to compute.
;       x: The domain of values to be spanned by the output polynomials.
; KEYWORD PARAMETERS:
;       domain: Explicity specify domain onto which to rescale the interval [-1,1]
; OUTPUTS:
;       Result: dblarr(nx, n+1).  The ith row contains the ith degree 
;               polynomial.
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       JEM, July, 2008. Written
;       JEM, Dec, 2010. Added domain keyword
;-

Function Cheb_array, n, x, domain=domain
  if (n eq 0) then return, x*0.+1
  if n_elements(domain) eq 0 then domain=minmax(x)
  min = domain[0]
  domainsize = domain[1]-min
  nx = n_elements(x)
  c = cheb_coeffs(n)
  xx = 2.*(x-min)/domainsize-1. ;rescale to [-1,1]
  result = dblarr(nx, n+1)
  for i = 0, n do result[*, i] = poly(xx, c[i, *])
  return, result
end
