FUNCTION minimize_me, P
  COMMON biweight_linfitblock, X, Y, b0, silent
  yfit = b0 + P[0] * X
  resid = biweight_location(Y - yfit, sigma)
  b0 += resid
  if ~keyword_set(silent) then print, sigma, b0, P[0]
  RETURN, sigma
END

FUNCTION biweight_linfit, Xin, Yin, sigma, weights, silent=silent1
  COMMON biweight_linfitblock
  silent=keyword_set(silent1)
  X = Xin
  Y = Yin
  result0 = linfit(X,Y,yfit=yfit)
  m0 = result0[1]
  b0 = result0[0]
  slope = amoeba(0.0001, function_name='minimize_me', p0=m0, scale=m0/10.)
  yfit = b0 + slope[0] * X
  resid = biweight_location(Y - yfit, sigma, weights)
  b0 += resid
  RETURN, [b0,slope[0]]
END
