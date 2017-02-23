Function covariance, data, datamean=datamean
  s=size(data, /dim)
  nparam = s[0]
  nsample = s[1]
  datamean = total(data, 2)/nsample
  covar = fltarr(nparam, nparam)
  for i=0, nparam-1 do begin
     for j=0, i do begin
        covarij = total((data[i,*]-datamean[i]) * (data[j,*]-datamean[j]))
        covarij /= nsample-1
        covar[i,j]=covarij
        covar[j,i]=covarij
     endfor
  endfor
  return, covar
end
