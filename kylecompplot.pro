Pro kylecompplot
  metals = [0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05]
  set_plot, 'ps'
  device, /encapsulated, /color, bits_per_pixel=8
  device, xsize=5., ysize=5., /inches, xoffset=0, yoffset=0
  device, filename='KyleMassComp.eps'
  plot, [0], /nodata, xrange=[0.2, 1.2], yrange=[46., 50.], $
        xtitle=textoidl('i_{775}-z_{850}'), $
        ytitle=textoidl('z_{850} for galaxy with 1 M_{sun} in stars'), $
        title='z=1.00 galaxy with star-formation burst'
  for zf = 2, 6 do begin
     izs = fltarr(6, 2)
     for imetal = 0, 5 do begin
        age = galage( 1., zf )
        origmass = 1./bc03starfrac( age, metals[imetal] )
        spec = bc03spec( age, metals[imetal], /erg, origmass = origmass )
        wave = spec->wavelength()
        lumdens = spec->flux()
        izs[imetal,*] = izmags4( wave, lumdens, 1.0 )
     endfor
     oplot, izs[*,0]-izs[*,1], izs[*,1], ps=symcat(16)
     oplot, izs[*,0]-izs[*,1], izs[*,1]
  endfor
  device, /close_file
  set_plot, 'X'
end
