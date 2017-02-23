Function Gini, X, abs=abs, err=err
  X=X[sort(X)]
  n=n_elements(X)
  if keyword_set(abs) then begin
     gini=1./(mean(abs(X))*n*(n-1)) * total( (2*findgen(n)-n-1)*abs(X) )
  endif else begin
     gini=1./(mean(X)*n*(n-1)) * total( (2*findgen(n)-n-1)*X )
  endelse
  if arg_present(err) then begin
     Gs=fltarr(200)
     for i=0,199 do begin
        w=floor(randomu(seed, n)*n)
        X1=X[w]
        Gs[i]=Gini(X1, abs=keyword_set(abs))
     endfor
     gini=biweight_location(Gs)
     err=biweight_scale(Gs-gini, /zero, /nan)
  endif
  return, gini
end
