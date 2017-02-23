Function Fit_RS, mags, colors, errs, fit=fit
  if n_elements(fit) eq 0 then fit=indgen(n_elements(mags))
  nobj=n_elements(mags)
  nfit=n_elements(fit)
  b=fltarr(5000)
  m=fltarr(5000)
  s=fltarr(5000)
  ;; so i'm doing the bootstrapping on just the "fit" sample
  ;; but i'll get the sigma from the entire sample later...
  for iter=0, 4999 do begin
;     w=floor(randomu(seed, nobj)*nobj)
;     mags1 = mags[w]  ;; sampling with replacement
;     colors1 = colors[w]
;     errs1 = errs[w]
     w=floor(randomu(seed, nfit)*nfit)
     mags1 = mags[fit[w]]  ;; sampling with replacement
     colors1 = colors[fit[w]]
     errs1 = errs[fit[w]]

     result1 = linfit( mags1, colors1, measure_errors=errs1, yfit=yfit )  ;fit
     b[iter] = result1[0]+result1[1]*23  ;;record results
     m[iter] = result1[1]
     resids = colors1-yfit
;     resids = colors-(result1[0]+result1[1]*mags)
     s1 = robust_sigma(resids, /zero)
     s[iter]=s1
  endfor
  bbar = mean(b)  ;; take mean of samples
  mbar = mean(m)
  sbar = mean(s)
  bsig = stdev(b)
  msig = stdev(m)
  ssig = stdev(s)
  yfit = bbar + mbar*(mags-23)
  resids = colors-yfit
  sig0 = robust_sigma(resids)  ;; sig0 is robust_sigma for average best-fit line
  errscat = errscat(errs)  ;; photometric expectation of scatter.
  @define_structs
  out = RSfit0
  out.slope = mbar
  out.slope_err = msig
  out.intercept = bbar
  out.intercept_err = bsig
  out.mscatter = sbar
  out.mscatter_err = ssig
  out.mscatter_best = sig0
  out.scatter_corr = errscat
  out.ngals = n_elements(mags)
  out.iscatter = sqrt(out.mscatter_best^2 - out.scatter_corr^2)
  return, out
end
