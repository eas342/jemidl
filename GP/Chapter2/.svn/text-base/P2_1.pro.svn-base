Function covarfn, x
  n=n_elements(x)
  covar=fltarr(n,n)
  for i=0, n-1 do begin
     for j=0, i do begin
        covar[i,j] = exp(-0.5*(x[i]-x[j])^2)
        covar[j,i] = covar[i,j]
     endfor
     covar[i,i] += 1.d-6 ;; for numerical stability
  endfor
  return, covar
end

;; covar is covariance of [x, xtest] vector
;; returns, ymean, ycovar
Function GPpredict, x, y, xtest, covar
  nx = n_elements(x)
  nxtest = n_elements(xtest)
  cxx = covar[0:nx-1, 0:nx-1]
  cxx_inv = invert(cxx)
  cxxtest = covar[0:nx-1, nx:*]
  cxtestxtest = covar[nx:*,nx:*]

  mean = transpose(cxxtest ## cxx_inv ## y)
  covar = cxtestxtest - cxxtest ## cxx_inv ## transpose(cxxtest)
  return, {mean:mean, covar:covar}
end

Pro P2_1
  common p2_1block, seed
  xrange=[-5,5]
  n=100
  x=dindgen(n)/(n-1)*(xrange[1]-xrange[0])+xrange[0]
  covar=covarfn(x)
  !p.multi=[0,1,2]
  plot, [0], /nodata, xrange=xrange, yrange=[-4,4]
  for i=0, 9 do begin
     y=mrandomn(seed,covar)
     oplot, x, y
  endfor

  data = {x:[-2,2,3], y:[-2,2,1]}
  covar = covarfn([data.x,x])
  post = GPpredict(data.x, data.y, x, covar)
  plot, data.x, data.y, xrange=xrange, yrange=[-4,4], ps=1, thick=2
  oplot, x, post.mean, thick=2
  for i=0, 9 do begin
     y=mrandomn(seed,post.covar)+post.mean
     oplot, x, y
  endfor
  !p.multi=0
end
