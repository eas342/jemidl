Function meanshiftkernel, x, kernelsize
  return, exp(-(x/kernelsize)^2)
end

Function meanshift, data, kernelsize, niter
  data1 = data
  for k=0, niter-1 do begin
     newdata = data1*0
     for i=0L, n_elements(data)-1 do begin
        newdata[i] = total(data*meanshiftkernel(data-data[i], kernelsize)) $
                     /total(meanshiftkernel(data-data[i], kernelsize))
     endfor
     data1 = newdata
  endfor
  return, newdata
end
