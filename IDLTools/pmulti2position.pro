Function PMulti2Position, region, position
  xrange = region[2]-region[0]
  xstart = region[0]
  yrange = region[3]-region[1]
  ystart = region[1]
  p = [xstart+position[0]*xrange, $
       ystart+position[1]*yrange, $
       xstart+position[2]*xrange, $
       ystart+position[3]*yrange]
  return, p
end
