Function Centroid, image, mask
  sz=size( image, /dim )
  make_2d, indgen(sz[0]), indgen(sz[1]), xx, yy
  xmoment = total(xx*image*mask)
  ymoment = total(yy*image*mask)
  bright = total(image*mask)
  xmoment /= bright
  ymoment /= bright
  return, [xmoment, ymoment]
end
