Function weightedmean, x, w, sigma=xsigma
  xbar = total(x*w)/total(w)
  n = n_elements(w)
  xsigma = 1./total(w) * 1./(n-1) * total((x-xbar)^2*w)
  return, xbar
end
