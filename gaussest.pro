Function gaussestloglikelihood, theta
  common gaussest_block, x, dx
  mu = theta[0]
  sigma = theta[1]
  if sigma  lt 0. then return, -1.0*!values.f_infinity

  sigmaext = sigma^2+dx^2
;  return, 1./sqrt(2*!pi*total(sigmaext))*exp(-0.5*total((mu-x)^2/sigmaext))
  return, -0.5*total(alog(2*!pi*sigmaext)) - total(0.5*(mu-x)^2/sigmaext)
end

;; Given data and uncertainties on that data,
;; find parent gaussian distribution with mean
;; and intrinsic dispersion that best fits the data.
Function gaussest, x1, dx1, covar=covar
  common gaussest_block
  x=x1
  dx=dx1
  mu_guess = mean(x)
  var_x = variance(x)
  sigma_guess = sqrt(var_x-median(dx)^2 > 0.1*var_x)
  mu_var = var_x/n_elements(x)
  sigma_var = var_x/(2*n_elements(x))
  parinfo = replicate( { name:'', $
                         jstargoal:20.0, $
                         rgoal:0.01, $
                         fixed:0 }, 2 )
  parinfo.name = ['mu','sigma']
;;  parinfo[0].fixed=1
  result = metropolis3( 'gaussestloglikelihood', [mu_guess, sigma_guess], $
                        diag_matrix([mu_var,sigma_var]), $
                        ntrial=1, trialniter=150L, parinfo=parinfo, $
                        miniter=100L, maxiter=100000L, nstop=4, /silent, $
                        paramchain=paramchain, valchain=valchain, fval=fval)
  bestval = max(valchain, wbest)
  wlikely=where(valchain gt bestval-alog(100.))
  covar = covariance(paramchain[*,wlikely])
  return, result
end
