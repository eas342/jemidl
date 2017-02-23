;;assume that abscissae are strictly monotonically increasing
Function interpextrap, ordinate, abscissa, var
  winterp = where(var ge min(abscissa) and var le max(abscissa))
  whigh = where(var gt max(abscissa))
  wlow = where(var lt min(abscissa))
  out = var*0
  if winterp[0] ne -1 then begin
     out[winterp] = interpol( ordinate, abscissa, var[winterp] )
  endif
  if whigh[0] ne -1 then begin
     n = n_elements(ordinate)
     slope = (ordinate[n-1]-ordinate[n-2])/(abscissa[n-1]-abscissa[n-2])
     dx = var[whigh]-abscissa[n-1]
     out[whigh] = ordinate[n-1]+slope*dx
  endif
  if wlow[0] ne -1 then begin
     n = n_elements(ordinate)
     slope = (ordinate[1]-ordinate[0])/(abscissa[1]-abscissa[0])
     dx = var[wlow]-abscissa[0]
     out[wlow] = ordinate[0]+slope*dx
  endif
  return, out
end
