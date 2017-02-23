;fix the slope, do grid search to minimize # RS scatter within 0.25 mag
;by 0.005 increments in intercept?

Function Fit_RS5, mags, colors, errs, fit=fit, fixslope=fixslope, int23=int23
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
  
  n = 0.4/0.005+1
  ints = findgen(n)/(n+1)+int23-0.2
  scat = ints*0.
  for i=0, n_elements(ints)-1 do begin
     int23try = ints[i]
     model = int23try+fixslope*(mags1-23)
     resid = colors1-model
     iscat = 0.03
     for j=0, 2 do begin
        wgal=where(abs(resid/sqrt(iscat^2+errs^2)) lt 2.5 and abs(resid) lt 0.1 or abs(resid) lt 0.06)
        if n_elementS(wgal) lt 5 then begin
           iscat=1.0 
           break
        endif
        measured_scatter = biweight_scale( resid[wgal], /zero )
        measurement_scatter = errscat2(errs[wgal])
        iscat = sqrt(measured_scatter^2 - measurement_scatter^2)
     endfor
     scat[i] = iscat
  endfor
  best = min( scat, /nan )
  wmin = where( scat eq best )
  int23out = mean(ints[wmin])

  model = int23out+fixslope*(mags-23)
  resid = colors-model
;  iqr = percentile(resid,[0.25,0.75])
;  select = where(resid gt iqr[0]-1.5*(iqr[1]-iqr[0]) and resid lt
;  iqr[1]+(iqr[1]-iqr[0]))
  select = where(resid gt -100.)

;  measured_scatter = robust_sigma( resid[select], /zero )
  measured_scatter = biweight_scale( resid[select], /zero )
  measurement_scatter = errscat2(errs[select])
  iscatter = sqrt(measured_scatter^2 - measurement_scatter^2)


  @define_structs
  out = RSfit0
  out.slope=fixslope
  out.slope_err=0.
  out.intercept=int23out
  out.intercept_err=0.
  out.mscatter=measured_scatter
  out.mscatter_err=0.
  out.mscatter_best=measured_scatter
  out.scatter_corr=measurement_scatter
  out.ngals=nobj
  out.iscatter=iscatter
  
  return, out
end
