;
;return a single asymmetry measurement w/o background asymmetry subtraction
;
Function Asym1, image, mask, center, interp=interp
  if keyword_set(interp) then begin
     rotimage = rot( image, 180.d, 1., $
                     center[0], center[1], /interp, /pivot )
  endif else begin
     rotimage = rot( image, 180.d, 1., $
                     center[0], center[1], cubic=-0.5, /pivot )
  endelse
;  rotresid = (image-rotimage)*mask
;  atv, rotresid
;  return, total(abs(rotresid))
  return, total(abs(image - rotimage)*mask)
end

;
;efficiently return a 3x3 array of asymmetry measurements by shifting 
;the center of rotation.
;
Function Asym9, image, mask, center, prevcenter, prev9, interp=interp
  if n_elements(prevcenter) eq 0 then begin
     asym9 = fltarr(3,3)
     for i=-1, 1 do begin
        for j=-1, 1 do begin
           asym9[i+1, j+1] = asym1( image, mask, center+[i,j]*0.2, $
                                    interp=interp )
        endfor
     endfor
     return, asym9
  endif
  asym9 = prev9
  case 1 of
     (center[0] lt prevcenter[0]) : begin 
        asym9 = shift( asym9, 1, 0 )
        asym9[0, *] = !values.f_nan
     end
     (center[0] gt prevcenter[0]) : begin
        asym9 = shift( asym9, -1, 0 )
        asym9[2, *] = !values.f_nan
     end
     else : 
  endcase
  case 1 of 
     (center[1] lt prevcenter[1]) : begin
        asym9 = shift( asym9, 0, 1 )
        asym9[*,0] = !values.f_nan
     end
     (center[1] gt prevcenter[1]) : begin
        asym9 = shift( asym9, 0, -1 )
        asym9[*,2] = !values.f_nan
     end
     else :
  endcase
  w=where(~finite(asym9))
  for iw=0, n_elements(w)-1 do begin
     ind = array_indices(asym9, w[iw])
     asym9[w[iw]] = asym1( image, mask, center+(ind-1)*0.2, $
                           interp=interp )
  endfor
  return, asym9
end

;
;Find a galaxy's asymmetry
;
Function Asymmetry, image, mask1, sky_sig, $
                    center=center, err=err, rotresid=rotresid, $
                    interp=interp, symmetric=symmetric
  interp=keyword_set(interp)
  if ~finite(total(image)) then begin
     err=!values.f_nan
     if arg_present(rotresid) then rotresid=image*!values.f_nan
     return, !values.f_nan
  endif
  s=size(image, /dim)
  if n_elements(center) eq 0 then center = s/2

  mask = mask1
  if keyword_set(symmetric) then begin
     mask = mask and ( rot( mask, 180.d, 1., $
                            center[0], center[1], $
                            /interp, /pivot ) gt 0.5 )
  endif

  asym9 = asym9( image, mask, center, $
                 interp=interp )
  num1 = min( asym9, wm )
  shift = array_indices( asym9, wm )

  iter = 0
  while total(shift eq [1,1]) ne 2 and iter lt 50 do begin
     prevcenter = center
     center += (shift-1)*0.2
     prev9 = asym9
     if keyword_set(symmetric) then begin
        mask = mask1 and ( rot( mask1, 180.d, 1., $
                                center[0], center[1], $
                                /interp, /pivot ) gt 0.5 )
     endif
     asym9 = asym9( image, mask, center, prevcenter, prev9, $
                    interp=interp )
     num1 = min( asym9, wm )
     shift = array_indices( asym9, wm )
     iter += 1
  endwhile

  if arg_present(rotresid) then begin
     if keyword_set(interp) then begin
        rotimage = rot( image, 180.d, 1., $
                        center[0], center[1], /interp, /pivot )
    endif else begin
        rotimage = rot( image, 180.d, 1., $
                        center[0], center[1], cubic=-0.5, /pivot )
     endelse
     rotresid = (image - rotimage)*mask
  endif
  num2arr = fltarr(21)
  for inum2=0, n_elements(num2arr)-1 do begin
     noise = randomn(seed, size(image,/dim))*sky_sig
     noise25 = fltarr(5,5)
     for i=-2,2 do begin
        for j=-2,2 do begin
           noise25[i+2,j+2] = asym1( noise, mask, center+[i,j]*0.2, $
                                     interp=interp )
        endfor
     endfor
     num2arr[inum2] = min(noise25)
  endfor

  num2 = median(num2arr)
  den = total(image*mask)
  asym=(num1-num2)/den
  err = biweight_scale(num2arr, /nan)/den
;  err = sqrt((moment(num2arr))[1])/den
;  print, string(format='(%"(%f - %f)/%f = %f")', num1, num2, den,
;                                 (num1-num2)/den)
  return, asym
end
