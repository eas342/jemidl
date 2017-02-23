Pro doMeasurePeakCenters, file
  if n_elements(file) ne 0 then $
     readcol, file, fileNames, shifts, format='A,F', /silent
  files=file_search('s??-*.fits.gz') ;;either std- or sci-
  for ifile=0, n_elements(files)-1 do begin
     if strpos(files[ifile], 'newwave') ne -1 then continue ;; ignore newwave files.
     if n_elements(file) ne 0 then begin
        w=where(files[ifile] eq fileNames)
        if w[0] ne -1 then begin
           measurePeakCenters, files[ifile], shifts[w[0]]
        endif else begin
           measurePeakCenters, files[ifile], 0.0
        endelse
     endif else begin
        measurePeakCenters, files[ifile], 0.0
     endelse
  endfor
end
