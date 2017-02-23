Function Fit_RS2, mags, colors, errs, fit=fit, iscatter=iscatter
  if n_elements(iscatter) eq 0 then iscatter=0.03
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
  result = linfit( mags1, colors1, $
                   measure_errors=sqrt(errs1^2+iscatter^2) )
  
  resid = colors - (result[0]+result[1]*mags)
  measured_scatter = stdev(resid)
  
  measurement_scatter = errscat2(errs)

  iscatter = sqrt( measured_scatter^2 - measurement_scatter^2 )
  @define_structs
  out = RSfit0
  out.slope = result[1]
  out.slope_err=0
  out.intercept=result[0]+result[1]*23
  out.intercept_err=0
  out.mscatter=measured_scatter
  out.mscatter_err=0
  out.mscatter_best=measured_scatter
  out.scatter_corr=measurement_scatter
  out.ngals=nobj
  out.iscatter=iscatter


  return, out
end
