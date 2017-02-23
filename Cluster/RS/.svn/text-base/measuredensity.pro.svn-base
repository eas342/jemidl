Function measuredensity, gals
  mags = findgen(18)/17*6+19.
  colors = findgen(20)/19.*2-0.5
  dens = fltarr( n_elements(mags)-1, n_elements(colors)-1 )
  for imag=0, n_elements(mags)-2 do begin
     for icolor=0, n_elements(colors)-2 do begin
        dens[imag,icolor] = total( gals.zmag gt mags[imag] $
                                   and gals.zmag le mags[imag+1] $
                                   and gals.iz gt colors[icolor] $
                                   and gals.iz le colors[icolor+1] )
     endfor
  endfor
  return, dens
end
