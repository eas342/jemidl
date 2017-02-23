Function Concentration, image, center, radius
  s = size(image,/dim)
  maxr = min(s)/2
  r = findgen(maxr)+1
  flux = apphot( image, center[0]+0.5, center[1]+0.5, r )
  fluxtot = interpol( flux, r, radius*1.5 )
  r80 = interpol( r, flux, fluxtot*0.8 )
  r20 = interpol( r, flux, fluxtot*0.2 )
  return, 5*alog10(r80/r20)
end
