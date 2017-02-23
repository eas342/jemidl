Function Mtot_M20, p
  common JEM$_M20, image, mask, xs, ys
  return, total(image*mask*((xs-p[0])*(xs-p[0])+(ys-p[1])*(ys-p[1])))
end

Function M20, image1, mask1, guess=guess
  common JEM$_M20
  image=image1
  mask=mask1
  s=size(image,/dim)
  xs=findgen(s[0])
  ys=findgen(s[1])
  make_2d, xs, ys

  if n_elements(guess) eq 0 then guess = s/2

  scale = [3,3]
  center = downhillsimplex( 'mtot_m20', guess, scale, tol=0.00001, fval=fval )

  ftot = total(image*mask)

  image=image[where(mask)]
  xs=xs[where(mask)]
  ys=ys[where(mask)]

  s=sort(image)
  fcum = total( (image)[reverse(s)], /cumulative )
  w=where( fcum lt 0.2*ftot)
  if w[0] eq -1 then return, !values.f_nan
  Mitot = total(image[w]*((xs[w]-center[0])^2+(ys[w]-center[1])^2))
  return, alog10(Mitot/fval[0])
end
