Function smoothness, image, skysigma, radius, center=center
  if radius lt 2. then return, !values.f_nan
  s=size( image, /dim )
  ksize = floor(min([s,radius]))>3
  xs=indgen(ksize)
  ys=indgen(ksize)
  sigma = radius*0.25/2.35482
  gauss = sampledgauss2d(xs, ys, [0., 1., sigma, sigma, ksize/2.-0.5, ksize/2.-0.5])
  gauss /= total(gauss)
  smooth = convol( image, gauss, /edge_truncate )
  resid = abs(image-smooth)
  dist_ellipse, dist, s, center[0], center[1], 1., 0.
  wnum = where(dist gt 0.25*radius and dist lt 1.5*radius)
  wden = where(dist lt 1.5*radius)
  if wnum[0] eq -1 or wden[0] eq -1 then return, !values.f_nan
  skyimage = randomn( seed, s[0], s[1] )*skysigma
  smoothsky = convol( skyimage, gauss, /edge_truncate )
  skyresid=abs(skyimage-smoothsky)
  num1=total(resid[wnum])
  num2=total(skyresid[wnum])
  den=total(image[wden])
  return, (num1-num2)/den
end
