; This is Josh Meyers' mod of amoeba.pro to implement
; minimization of an object method instead of a function.

; $Id: amoeba.pro,v 1.14 2005/02/01 20:24:12 scottm Exp $
; Copyright (c) 1996-2005, Research Systems, Inc.  All rights reserved.
;	Unauthorized reproduction prohibited.


Function amotryobj, p, y, psum, func, obj, ihi, fac
; Extrapolates by a factor fac through the face of the simplex, across
; from the high point, tries it and replaces the high point if the new
; point is better.

  compile_opt hidden
  fac1 = (1.0 - fac) / n_elements(psum)
  fac2 = fac1  - fac
  ptry = psum * fac1 - p[*,ihi] * fac2
  ytry = call_method(func, obj, ptry)  ;Eval fcn at trial point
  if ytry lt y[ihi] then begin    ;If its better than highest, replace highest
    y[ihi] = ytry
    psum = psum + ptry - p[*,ihi]
    p[0,ihi] = ptry
    endif
  return, ytry
end


Function Amoebaobj, ftol, FUNCTION_NAME=func, OBJECT=obj, $
                    FUNCTION_VALUE=y, NCALLS = ncalls, NMAX = nmax, $
                    P0 = p0, SCALE=scale, SIMPLEX=p
; The Numerical Recipes procedure Amoeba, with some embellishments.
;
; Reference: Numerical Recipes, 2nd Edition, Page 411.
;  P = fltarr(ndim, ndim+1)
;+
; NAME:
;	AMOEBAOBJ
;
; PURPOSE:
;	Multidimensional minimization of a method FUNC(X), where
;	X is an N-dimensional vector, using the downhill simplex
;	method of Nelder and Mead, 1965, Computer Journal, Vol 7, pp 308-313.
;
;	This routine is based on the AMOEBA routine, Numerical
;	Recipes in C: The Art of Scientific Computing (Second Edition), Page
;	411, and is used by permission.
;
; CATEGORY:
;	Function minimization/maximization. Simplex method.
;
; CALLING SEQUENCE:
;	Result = AMOEBAOBJ(Ftol, ....)
; INPUTS:
;    FTOL:  the fractional tolerance to be achieved in the method
;	value.  e.g. the fractional decrease in the method value in the
;	terminating step.  This should never be less than the
;	machine's single or double precision.
; KEYWORD PARAMETERS:
;    FUNCTION_NAME: a string containing the name of the method to
;	be minimized.  If omitted, the method FUNC is minimized.
;	This method must accept an Ndim vector as its only parameter and
;	return a scalar single or double precision floating point value as its
;	result. 
;    FUNCTION_VALUE: (output) on exit, an Ndim+1 element vector
;	containing the method values at the simplex points.  The first
;	element contains the method minimum. 
;    NCALLS: (output) the of times the method was evaluated. 
;    NMAX: the maximum number of method evaluations allowed
;	before terminating.  Default = 5000.
;    P0: Initial starting point, an Ndim element vector.  The starting
;	point must be specified using either the keyword SIMPLEX, or P0 and
;	SCALE.  P0 may be either single or double precision floating.
;	For example, in a 3-dimensional problem, if the initial guess
;	is the point [0,0,0], and it is known that the method's
;	minimum value occurs in the interval: -10 <
;	X(0) < 10, -100 < X(1) < 100, -200 < X(2) < 200, specify: P0=[0,0,0],
;	SCALE=[10, 100, 200]. 
;    SCALE: a scalar or Ndim element vector contaiing the problem's
;	characteristic length scale for each dimension.
;	SCALE is used with P0 to form an initial (Ndim+1) point simplex.
;	If all dimensions have the same	scale, specify a scalar.
;    SIMPLEX: (output and/or optional input) On input, if P0 and SCALE
;	are not set, SIMPLEX contains the Ndim+1 vertices, each of
;	Ndim elements, of starting simplex, in either single or double
;	precision floating point, in an (Ndim, Ndim+1) array. On output,
;	SIMPLEX contains the simplex, of dimensions (Ndim, Ndim+1), enclosing
;	the method minimum.  The first point, Simplex(*,0), corresponds to
;	the method's minimum.
;
; OUTPUTS:
;   Result: If the minimum is found, an Ndim vector, corresponding to
;	the Method's minimum value is returned.  If a method minimum
;	within the given tolerance, is NOT found in the given number of
;	evaluations, a scalar value of -1 is returned.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	None.
;
; PROCEDURE:
;	This procedure implements the Simplex method, described in
;	Numerical Recipes, Section 10.4.  See also the POWELL procedure.
;
;	Advantages:  requires only method evaluations, not
;	derivatives, may be more reliable than the POWELL method.
;	Disadvantages: not as efficient as Powell's method, and usually
;	requires more method evaluations.
;
;	Results are performed in the mode (single or double precision)
;	returned by the user-supplied method.  The mode of the inputs P0,
;	SCALE, or SIMPLEX, should match that returned by the method. The
;	mode of the input vector supplied to the user-written method, is
;	determined by P0, SCALE, or SIMPLEX.
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;	DMS, May, 1996.	Written.
;       JEM, Jul, 2008. Changed to use object methods instead of functions...
;-

if keyword_set(scale) then begin  ;If set, then p0 is initial starting pnt
   ndim = n_elements(p0)
   p = p0 # replicate(1.0, ndim+1)
   for i=0, ndim-1 do p[i,i+1] = p0[i] + scale[i < (n_elements(scale)-1)]
endif

s = size(p)
if s[0] ne 2 then message, 'Either (SCALE,P0) or SIMPLEX must be initialized'
ndim = s[1]			;Dimensionality of simplex
mpts = ndim+1			;# of points in simplex
if n_elements(func) eq 0 then func = 'FUNC'
if n_elements(nmax) eq 0 then nmax = 5000L

y = replicate(call_method(func, obj, p[*,0]), mpts)  ;Init Y to proper type
for i=1, ndim do y[i] = call_method(func, obj, p[*,i]) ;Fill in rest of the vals
ncalls = 0L
psum = total(p,2)

while ncalls le nmax do begin		;Each iteration
  s = sort(y)
  ilo = s[0]		;Lowest point
  ihi = s[ndim]		;Highest point
  inhi = s[ndim-1]	;Next highest point
  d = abs(y[ihi]) + abs(y[ilo]) ;Denominator = interval
  if d ne 0.0 then rtol = 2.0 * abs(y[ihi]-y[ilo])/d $
  else rtol = ftol / 2.         ;Terminate if interval is 0
  
  if rtol lt ftol then begin ;Done?
      t = y[0] & y[0] = y[ilo] & y[ilo] = t ;Sort so fcn min is 0th elem
      t = p[*,ilo] & p[*,ilo] = p[*,0] & p[*,0] = t
      return, t                 ;params for fcn min
  endif
    
  ncalls = ncalls + 2
  ytry = amotryobj(p, y, psum, func, obj, ihi, -1.0)
  if ytry le y[ilo] then ytry = amotryobj(p,y,psum, func, obj, ihi, 2.0) $
  else if ytry ge y[inhi] then begin
    ysave = y[ihi]
    ytry = amotryobj(p,y,psum,func, obj, ihi, 0.5)
    if ytry ge ysave then begin
	for i=0, ndim do if i ne ilo then begin
	  psum = 0.5 * (p[*,i] + p[*,ilo])
	  p[*,i] = psum
	  y[i] = call_method(func, obj, psum)
	  endif
	ncalls = ncalls + ndim
	psum = total(p,2)
	endif		;ytry ge ysave
    endif else ncalls = ncalls  - 1
  endwhile
  return, -1		;Here, the function failed to converge.
end
