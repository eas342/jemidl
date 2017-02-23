Function wave2bins, wave
  binbounds = 0.5*(wave[1:*]+wave[0:n_elements(wave)-2])
  binbounds = [2*wave[0] - binbounds[0], $
               binbounds, $
               2*wave[n_elements(wave)-1]-binbounds[n_elements(binbounds)-1]]
  return, binbounds
end
