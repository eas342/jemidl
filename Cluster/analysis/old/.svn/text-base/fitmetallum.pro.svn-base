Function fitmetallum, specgals, ap_radius=ap_radius, ap_factor=ap_factor, _extra=extra
  if n_elements(ap_factor) eq 0 and n_elements(ap_radius) eq 0 then ap_factor=1.

  diameters=[1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0, $
             22.0,24.0,26.0,28.0,30.0,32.0,40.0,60.0,80.0,120.0,160.0,220.0]

;  wno2=where( strpos(specgals.comment, 'oII') eq -1 $
;              and abs(specgals.velocity) lt 2*(specgals.vel_disp) )

  for i=0, n_elements(specgals)-1 do begin
     if n_elements(ap_factor) ne 0 then begin
        imag = interpol( specgals[i].imag_aper, diameters, 2.*ap_factor*specgals[i].zr_e )
        zmag = interpol( specgals[i].zmag_aper, diameters, 2.*ap_factor*specgals[i].zr_e )
     endif else begin
        imag = interpol( specgals[i].imag_aper, diameters, 2.*ap_radius )
        zmag = interpol( specgals[i].zmag_aper, diameters, 2.*ap_radius )
     endelse
     zmagbest = specgals[i].zmag_best
     imag = zmagbest + imag-zmag
     zmag = zmagbest

     metallum1 = fitgal( imag, zmag, specgals[i].z, zcluster=specgals[i].zcluster, _extra=extra )
     print, metallum1.metal, metallum1.lum, format='(2E13.4)'
     if n_elements(metallum) eq 0 then metallum=metallum1 $
     else metallum=[metallum, metallum1]
  endfor
  return, metallum
end
