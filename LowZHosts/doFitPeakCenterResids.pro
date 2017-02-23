Pro doFitPeakCenterResids
  files=file_search('wavecal/*.peakCenters.dat')
  for ifile=0, n_elements(files)-1 do begin
     fitPeakCenterResids, files[ifile]
  endfor
end
