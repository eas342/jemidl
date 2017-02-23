Pro doFitSkyWarp
  files=file_search('sci-?1108??_????.fits.gz')
  for ifile=0, n_elements(files)-1 do begin
     shift=0.0
     if files[ifile] eq 'sci-r110830_0044.fits.gz' then shift=86.2
     if files[ifile] eq 'sci-r110830_0045.fits.gz' then shift=86.2
     if files[ifile] eq 'sci-r110830_0046.fits.gz' then shift=95.3
     fitSkyWarp, files[ifile], shift
  endfor
end
