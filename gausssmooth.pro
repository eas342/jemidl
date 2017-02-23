Function gausssmooth, flux, sigma
  n = fix(sigma*3)*2+1
  x = findgen(n)
  y = exp(-(x-(n-1)/2)^2/(2*sigma^2))
  y /= total(y)
  return, convol( flux, y, /center, /edge_trunc )
end
