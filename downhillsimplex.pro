Function getval_downhillsimplex, fcn, obj, abscissa
  if ~obj_valid(obj) then begin
     return, call_function( fcn, abscissa )
  endif else begin
     return, call_method( fcn, obj, abscissa )
  endelse
end

Pro improve, fcn, simplex, fval, object=object
  ndim = n_elements(fval)-1
  ;; get midpoint of best 'side'
  mid = total(simplex[0:ndim-1,*],1)/ndim
  ;; get reflection point
  reflect = mid + (mid - simplex[ndim,*])
  ;; get extension point
  extend = reflect + (reflect - mid)
  ;; get fcn values
  fE = getval_downhillsimplex( fcn, object, extend )
  fW = fval[ndim]
;  fW = getval_downhillsimplex( fcn, object, simplex[ndim,*] )  ;this is potentially redundant
  ;; decision tree...
  ;; take extension if better
  if fE lt fW then begin
     simplex[ndim,*] = extend
     fval[ndim] = fE
     return
  endif
  ;; take reflection if better
  fR = getval_downhillsimplex( fcn, object, reflect )
  if fR lt fW then begin
     simplex[ndim,*] = reflect
     fval[ndim] = fR
     return
  endif
  ;; more trial abscissae. take them if better.
  c1 = 0.5 * (mid + simplex[ndim,*])
  fC1 = getval_downhillsimplex( fcn, object, c1 )
  if fC1 lt fW then begin
     simplex[ndim,*] = c1
     fval[ndim] = fC1
  endif
  c2 = 0.5 * (mid + reflect)
  fC2 = getval_downhillsimplex( fcn, object, c2 )
  if fC2 lt fW then begin
     simplex[ndim,*] = c2
     fval[ndim] = fC2
  endif
  ;; didn't find a better abscissa.  shrick simplex toward best abscissa
  for i=1, ndim do begin
     simplex[i,*] = 0.5 * (simplex[0,*] + simplex[i,*])
     fval[i] = getval_downhillsimplex( fcn, object, simplex[i,*] )
  endfor
end

Pro sortsimplex, fval, simplex
  wsort = sort(fval)
  fval = fval[wsort]
  tmp_simplex = simplex
  for i=0, n_elements(fval)-1 do tmp_simplex[i,*] = simplex[wsort[i],*]
  simplex = tmp_simplex
end

Function downhillsimplex, fcn, start, scale, $
                          tol=tol, object=object, $
                          miniter=miniter, maxiter=maxiter, fval=fval, verbose=verbose
  ;; initial setup
  ndim = n_elements(start)
  simplex = fltarr(ndim+1, ndim)
  fval = fltarr(ndim+1)
  if n_elements(object) eq 0 then object = ''
  if n_elements(maxiter) eq 0 then maxiter = 5000
  if n_elements(tol) eq 0 then tol = 0.01
  if n_elements(miniter) eq 0 then miniter = 1

  ;; fill in initial simplex and function values
  simplex[0,*] = start
  fval[0] = getval_downhillsimplex( fcn, object, simplex[0,*] )
  for i=0, ndim-1 do begin
     scale1 = fltarr(ndim)
     scale1[i] = scale[i]
     simplex[i+1,*] = start+scale1
  endfor
  for i=1, ndim do fval[i] = getval_downhillsimplex( fcn, object, simplex[i,*] )

  ;; sort simplex
  sortsimplex, fval, simplex
  
  for i=0, maxiter-1 do begin
     improve, fcn, simplex, fval, object=object
     sortsimplex, fval, simplex
     if i lt miniter then continue
     test = 2.0*(fval[ndim]-fval[0])/(fval[ndim]+fval[0])
     if keyword_set(verbose) then print, test
     if test lt tol then break
  endfor
  return, simplex[0,*]
end
