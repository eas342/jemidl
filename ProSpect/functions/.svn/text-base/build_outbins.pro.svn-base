Function Build_Outbins, $
   outbins=outbins, $
   dlam=dlam, $
   dloglam=dloglam, $
   wavemin=wavemin, $
   wavemax=wavemax

  if n_elements(outbins) ne 0 then return, outbins

  if n_elements(dlam) eq 0 and $
     n_elements(dloglam) eq 0 and $
     n_elements(outbins) eq 0 then begin
     message, 'output binning not specified.'
  endif
  if n_elements(dlam) ne 0 then begin
     range = wavemax - wavemin
     nbins = (range/dlam)+1
     outbins = dindgen(nbins)*dlam+wavemin
  endif else if n_elements(dloglam) ne 0 then begin
     range = alog10(1.0d * wavemax/wavemin)
     nbins = (range/dloglam)+1
     loutbins = dindgen(nbins)*dloglam+alog10(wavemin)
     outbins = 10.^loutbins
  endif
  return, outbins
end
