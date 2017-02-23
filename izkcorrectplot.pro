Pro izkcorrectplot, z, overplot=overplot, color=color
  ages = [1.2, 1.6, 2.0]*1.e9
  metals = [0.1, 0.5, 1.0, 2.0, 3.0]*0.0177


  if ~keyword_set(overplot) then plot, [0], /nodata, xrange=[0.,1.3], yrange=[0.,1.3]
  if z gt 1.237 then agediff = galage( 1.237, z ) else $
     agediff = -galage( z, 1.237 )
  print, ages-agediff
  print, agediff
  for iage=0, n_elements(ages)-1 do begin
     for imetal=0, n_elements(metals)-1 do begin
        obsspec = bc03spec( ages[iage]-agediff, metals[imetal], z=z )
        kspec = bc03spec( ages[iage], metals[imetal], z=1.237 )
        obsiz = izmags(obsspec, frame='obs')
        kiz = izmags(kspec, frame='obs')
        if n_elements(color) eq 0 then $
           oplot, [obsiz[0]-obsiz[1]], [kiz[0]-kiz[1]], ps=1 $
        else $
           oplot, [obsiz[0]-obsiz[1]], [kiz[0]-kiz[1]], ps=1, color=color
     endfor
  endfor
end
