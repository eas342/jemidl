Pro doApplyWaveCal
  files=file_search('s*.fits.gz')
  for ifile=0, n_elements(files)-1 do begin
     if strpos(files[ifile], 'newwave') ne -1 then continue
     print, files[ifile]
     hdr0=headfits(files[ifile],exten=0,/silent)
     targName = strtrim(sxpar(hdr0, 'TARGNAME'), 2)
     if strmid(targName,0,2) eq 'HD' then continue
     targName = lowzhost_nametranslate(targName)
     if sxpar(hdr0,'ELAPTIME') lt 10 then continue
     expName = strmid((strsplit(files[ifile],'.',/extract))[0],4)
     calFileName = 'wavecal/'+targName+'.'+expName+'.peakCenterResidFit.dat'
     applyWaveCal, files[ifile], calFileName
  endfor
end
