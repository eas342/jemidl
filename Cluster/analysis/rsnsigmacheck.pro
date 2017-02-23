Function rsnsigmacheck, $
   zmag, iz, iz_err, nsigma, magrange, bluelimit, $
   iscatter=iscatter, intercept=intercept, slope=slope

  if n_elements(intercept) eq 0 then intercept = 0.95
  if n_elements(slope) eq 0 then slope = -0.05
  if n_elements(iscatter) eq 0 then iscatter=0.05

  line = intercept + slope*(zmag-23)
  resid = iz-line
  nsig = resid/sqrt(iz_err^2+iscatter^2)
  return, ( nsig lt nsigma ) $
          and ( nsig gt -nsigma ) $
          and ( zmag gt magrange[0] ) $
          and ( zmag lt magrange[1] ) $
          and iz ge bluelimit
end
