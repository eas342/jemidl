Function Percentile, array, percentile
  fsort = sort(array)
  nsort = total(finite(array))
  return, array[fsort[floor(percentile*nsort)]]
end
