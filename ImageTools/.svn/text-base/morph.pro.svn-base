Function Morph, $
   image, $
   mask, $
   sky_sigma, $
   symmetric=symmetric, $
   rotresid=rotresid, $
   new_mask=new_mask, $
   _extra=extra
;possible extra: eta

  qpflux=quasipetroflux( image[where(mask)], _extra=extra )
  new_mask = mask and (image gt qpflux)
  if total(new_mask) le 2 then begin
     rotresid=image*!values.f_nan
     return, { qpflux:qpflux, gini:!values.f_nan, ginierr:!values.f_nan, $
               asym:!values.f_nan, asymerr:!values.f_nan }
  endif
  gini = gini( image[where(new_mask)], /abs, err=ginierr )
  if n_elements(center) eq 0 then center = centroid( image, mask )
  if arg_present(rotresid) then $
     asym = asymmetry( image, new_mask, sky_sigma, center=center, err=asymerr, $
                       symmetric=keyword_set(symmetric), rotresid=rotresid ) $
  else $
     asym = asymmetry( image, new_mask, sky_sigma, center=center, err=asymerr, $
                       symmetric=keyword_set(symmetric) )
  snr = mean(image[where(new_mask)]/sky_sigma)
  return, { qpflux:qpflux, gini:gini, ginierr:ginierr, asym:asym, asymerr:asymerr, snr:snr }
end
