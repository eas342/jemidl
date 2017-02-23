Function Fit_RS4_Posterior, slope
  return, Fit_RS4_Likelihood(slope)*Fit_RS4_prior(slope)
end

Function Fit_RS4_Likelihood, slope
  common JEM$_Fit_RS4_Likelihood, pslope, pwid, mags1, colors1, errs1
  resid = colors1-slope*mags1
  correction = total(resid/(errs1^2+0.03^2))/total(1./(errs1^2+0.03^2))
  resid -= correction
  return, total(-1.*resid^2/(2*(errs1^2+0.03^2)))
end

Function Fit_RS4_prior, slope
  common JEM$_Fit_RS4_Likelihood
  return, exp(-1.*(slope-pslope)^2/(2*(pwid^2)))
end

Function Fit_RS4, mags, colors, errs, fit=fit, $
                  pslope=pslope1, pwid=pwid1
  common JEM$_Fit_RS4_Likelihood
  pslope=pslope1
  pwid=pwid1

  w=where(finite(errs))
  mags=mags[w]
  colors=colors[w]
  errs=errs[w]

  nobj = n_elements(mags)  
  if n_elements(fit) eq 0 then fit = indgen(nobj)
  nfit = n_elements(fit)

  mags1 = mags[fit]
  colors1 = colors[fit]
  errs1 = errs[fit]
  
  slopes = findgen(50)/49.*(0.01+0.07)-0.07
  likes=slopes*0
  priors=slopes*0
  for i=0, 49 do begin
     likes[i] = fit_rs4_likelihood(slopes[i])
     priors[i] = fit_rs4_prior(slopes[i])
  endfor
  junk=fsc_color(/all, color=ctable)
  likes=exp(likes)
  plot, slopes, likes/max(likes)
  oplot, slopes, priors/max(priors), color=ctable.red
  oplot, slopes, priors*likes/(max(priors*likes)), color=ctable.blue
  legend, ['prior','like','post'], textcolor=[ctable.red,ctable.white, ctable.blue], $
          charsize=2, linestyle=[0,0,0], color=[ctable.red,ctable.white, ctable.blue]
  

  junk = max(priors*likes, m)
  slope = slopes[m]

  int23guess = 1.0
  resid = colors1 - (int23guess+slope*(mags1-23))
  correction = total(resid/errs1^2)/total(1./errs1^2)

  resid = colors - (1.+correction+slope*(mags1-23))
  measured_scatter = stdev(resid)
  
  measurement_scatter = errscat2(errs)

  iscatter = sqrt( measured_scatter^2 - measurement_scatter^2 )
  @define_structs
  out = RSfit0
  out.slope = slope
  out.slope_err=0
  out.intercept=1.+correction
  out.intercept_err=0
  out.mscatter=measured_scatter
  out.mscatter_err=0
  out.mscatter_best=measured_scatter
  out.scatter_corr=measurement_scatter
  out.ngals=nobj
  out.iscatter=iscatter

  return, out
end
