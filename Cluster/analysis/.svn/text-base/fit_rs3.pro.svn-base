;;fix the slope

Function Fit_RS3, mags, colors, errs, fit=fit, slope=slope
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

  int23guess = 1.0
  resid = colors1 - (int23guess+slope*(mags1-23))
  correction = total(resid/errs1^2)/total(1./errs1^2)

  resid = colors - (1.+correction+slope*(mags1-23))

  iqr = percentile(resid, [0.25, 0.75])
  select = where(resid gt iqr[0] - 1.5*(iqr[1]-iqr[0]) $
                 and resid lt iqr[1] + 1.5*(iqr[1]-iqr[0]))

;  measured_scatter = robust_sigma(resid[select], /zero)
  measured_scatter = biweight_scale(resid[select], /zero)
  measurement_scatter = errscat2(errs[select])
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
  out.ngals=n_elements(select)
  out.iscatter=iscatter

  return, out
end
