Function ConcRank, X
  n=n_elements(X)
  if n le 3 then return, !values.f_nan
  s=sort(X)
  t=total(X)
  cum = total(X[reverse(s)], /cumulative)
  n80 = interpol( indgen(n), X[s], 0.8*t )
  n20 = interpol( indgen(n), X[s], 0.2*t )
  return, 2.5*alog(n80/n20)
end
