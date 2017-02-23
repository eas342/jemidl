Function ensquaredenergySirianni, filter, side
  oversamplefactor = 9
  overside = oversamplefactor*201
  x = (findgen( overside ) - overside/2)/oversamplefactor
  y = (findgen( overside ) - overside/2)/oversamplefactor
  make_2d, x, y
  r = sqrt(x*x+y*y)
  array = encircledenergySirianni( filter, r, /diff )
  array /= total(array)
  w=where(abs(x) lt side/2 and abs(y) lt side/2)
  return, total(array[w])
end
