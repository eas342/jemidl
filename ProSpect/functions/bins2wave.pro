Function bins2wave, bins
  return, 0.5*(bins[1:*]+bins[0:n_elements(bins)-2])
end
