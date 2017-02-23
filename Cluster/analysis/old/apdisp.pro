Pro DoApDispPlots, specgals, clustername, ap_radius, ap_factor
  red = getcolor('red', !d.table_size-2)
  green = getcolor('green', !d.table_size-3)
  blue = getcolor('blue', !d.table_size-4)

  diameters=[1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0, $
             22.0,24.0,26.0,28.0,30.0,32.0,40.0,60.0,80.0,120.0,160.0,220.0]
  aps=[0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 13.0, 16.0, 20.0, 25., 30., 35., 40.]
  facs=[0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0]
  
  wno2 = where( strpos( specgals.comment, 'oII' ) eq -1 $
                and abs(specgals.velocity) lt 2.*specgals.vel_disp $
                and specgals.clustername eq clustername )
  
  ;;;;;;;;;;;;
  ; fixed ap

  sigmas = dblarr( n_elements(aps) )
  for iap=0, n_elements(aps)-1 do begin
     ap=2.*aps[iap]
     delvarx, mags
     for i=0, n_elements(wno2)-1 do begin
        imag = double(interpol( specgals[wno2[i]].imag_aper, diameters, ap ))
        zmag = double(interpol( specgals[wno2[i]].zmag_aper, diameters, ap ))
        zmagbest = specgals[wno2[i]].zmag_best
        imag = zmagbest + imag-zmag
        zmag = zmagbest
        if n_elements(mags) eq 0 then mags = {imag:imag, zmag:zmag} $
        else mags= [mags, {imag:imag, zmag:zmag}]
     endfor
     w = where( mags.zmag gt 19 $
                and mags.zmag lt 25 $
                and (mags.imag - mags.zmag) gt 0.5 $
                and (mags.imag - mags.zmag) lt 1.5 )
     fit = biweight_linfit( mags[w].zmag, (mags.imag - mags.zmag)[w], sigma )
     sigmas[iap]=sigma
  endfor
  minsig=min(sigmas,m)

  plot, aps, sigmas
  oplot, [ap_radius, ap_radius], [-1.e9, 1.e9], color=green
  xyouts, 15, max(sigmas)*0.8, string(minsig, format='(F7.5)')
  xyouts, 15, max(sigmas)*0.6, string(aps[m], format='(F4.1)')
  xyouts, 5, max(sigmas)*0.9, clustername
  

  ;;;;;;;;;;;;;;;;
  ; radius factor

  sigmas = dblarr( n_elements(facs) )
  for ifac=0, n_elements(facs)-1 do begin
     delvarx, mags
     for i=0, n_elements(wno2)-1 do begin
        ap=facs[ifac]*2.*specgals[wno2[i]].zr_e > 3.
        imag = double(interpol( specgals[wno2[i]].imag_aper, diameters, ap ))
        zmag = double(interpol( specgals[wno2[i]].zmag_aper, diameters, ap ))
        zmagbest = specgals[wno2[i]].zmag_best
        imag = zmagbest + imag-zmag
        zmag = zmagbest
        if n_elements(mags) eq 0 then mags = {imag:imag, zmag:zmag} $
        else mags= [mags, {imag:imag, zmag:zmag}]
     endfor
     w = where(mags.zmag gt 19 $
               and mags.zmag lt 25 $
               and (mags.imag - mags.zmag) gt -1 $
               and (mags.imag - mags.zmag) lt 2.)
     fit = biweight_linfit( mags[w].zmag, (mags.imag - mags.zmag)[w], sigma )     
     sigmas[ifac]=sigma
  endfor
  plot, facs, sigmas
  oplot, [ap_factor, ap_factor], [-1.e9, 1.e9], color=green
  xyouts, 2.5, max(sigmas)*0.8, string(min(sigmas,m), format='(F7.5)')
  xyouts, 2.5, max(sigmas)*0.6, string(facs[m], format='(F4.1)')
  
  ;;;;;;;;;;;;;
  ; CMD plot 1

  ap=ap_radius
  delvarx, mags
  for i=0, n_elements(wno2)-1 do begin
     imag = double(interpol( specgals[wno2[i]].imag_aper, diameters, ap ))
     zmag = double(interpol( specgals[wno2[i]].zmag_aper, diameters, ap ))
     zmagbest = specgals[wno2[i]].zmag_best
     imag = zmagbest + imag-zmag
     zmag = zmagbest
     if n_elements(mags) eq 0 then mags = {imag:imag, zmag:zmag} $
     else mags= [mags, {imag:imag, zmag:zmag}]
  endfor
  w = where( mags.zmag gt 19. $
             and mags.zmag lt 25. $
             and (mags.imag - mags.zmag) gt 0.5 $
             and (mags.imag - mags.zmag) lt 1.5 )
  fit = biweight_linfit( mags[w].zmag, (mags.imag - mags.zmag)[w], sigma, weights )

  resid = (mags.imag -mags.zmag)[w] - (fit[0]+fit[1]*mags[w].zmag)
  print, clustername, sqrt((moment(resid))[1])

  plot, mags[w].zmag, (mags.imag-mags.zmag)[w], xrange=[19,25], yrange=[0.5,1.5], ps=2, xstyle=1, ystyle=1
  w2 = where( weights ne 0. )
  oplot, mags[w[w2]].zmag, (mags.imag-mags.zmag)[w[w2]], ps=2, color=red
  xends = minmax(mags[w[w2]].zmag)
  oplot, xends, xends*fit[1]+fit[0], color=blue
  xyouts, 19.5, 1.32, string( fit[1], format='(F8.4)' )
  xyouts, 19.5, 1.2, string( fit[0]+21*fit[1], format='(F8.4)' )
  xyouts, 19.5, 1.08, string( sigma, format='(F7.5)' )


  ;;;;;;;;;;;;;
  ; CMD plot 2

  delvarx, mags
  for i=0, n_elements(wno2)-1 do begin
     ap=2.*ap_factor*specgals[wno2[i]].zr_e
     imag = double(interpol( specgals[wno2[i]].imag_aper, diameters, ap ))
     zmag = double(interpol( specgals[wno2[i]].zmag_aper, diameters, ap ))
     zmagbest = specgals[wno2[i]].zmag_best
     imag = zmagbest + imag-zmag
     zmag = zmagbest
     if n_elements(mags) eq 0 then mags = {imag:imag, zmag:zmag} $
     else mags= [mags, {imag:imag, zmag:zmag}]
  endfor
  w = where( mags.zmag gt 19. $
             and mags.zmag lt 25. $
             and (mags.imag - mags.zmag) gt 0.5 $
             and (mags.imag - mags.zmag) lt 1.5 )
  fit = biweight_linfit( mags[w].zmag, (mags.imag - mags.zmag)[w], sigma, weights )
  plot, mags[w].zmag, (mags.imag-mags.zmag)[w], xrange=[19,25], yrange=[0.5,1.5], ps=2, xstyle=1, ystyle=1
  w2 = where( weights ne 0. )
  oplot, mags[w[w2]].zmag, (mags.imag-mags.zmag)[w[w2]], ps=2, color=red
  xends = minmax(mags[w[w2]].zmag)
  oplot, xends, xends*fit[1]+fit[0], color=blue
  xyouts, 19.5, 1.32, string( fit[1], format='(F8.4)' )
  xyouts, 19.5, 1.2, string( fit[0]+22.*fit[1], format='(F8.4)' )
  xyouts, 19.5, 1.08, string( sigma, format='(F7.5)' )
end

Pro ApDisp, objs, ap_radius, ap_factor
  if n_elements(ap_radius) eq 0 then ap_radius=10.
  if n_elements(ap_factor) eq 0 then ap_factor=0.75
  colors = Obj_New("IDLgrPalette")
  colors->LoadCT, 2
  colors->GetProperty, Red=r, Green=g, Blue=b
  Obj_Destroy, colors
  TVLCT, r, g, b

  thisDevice=!d.name
  set_plot, 'ps'
  device, /encapsulated, color=1, bits_per_pixel=8, landscape=0
  device, xsize=7., ysize=9., /inches
  device, filename='apdisp.eps'
  
  specgals=specgals(objs)
  !p.multi = [0,4,6]
  DoApDispPlots, specgals, 'C', ap_radius, ap_factor
  DoApDispPlots, specgals, 'E', ap_radius, ap_factor
  DoApDispPlots, specgals, 'U', ap_radius, ap_factor
  DoApDispPlots, specgals, 'X', ap_radius, ap_factor
  DoApDispPlots, specgals, 'R', ap_radius, ap_factor
  DoApDispPlots, specgals, 'Y', ap_radius, ap_factor

  LoadCt, 0
  !p.multi=0
  device, /close_file
  set_plot, thisDevice
end
