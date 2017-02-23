Pro CRchecker, snname, setup, i_indiv
  dir='~/scp4/LowZHosts/JMreduction/'
  setups = ['B460','B560','R5000','R7500','R1000B','R2000B','R2500R']
  longsetups = ['blue600_4000_d460', $
                'blue600_4000_d560', $
                'red600_5000', $
                'red600_7500', $
                'R1000B', $
                'R2000B', $
                'R2500R']
  w=(where(setup eq setups))[0]
  checkfile=file_search(dir+'*/'+longsetups[w]+'/FinalCR/'+snname+'_check.fits')
  fluxfile=file_search(dir+'*/'+longsetups[w]+'/FinalCR/'+snname+'_F.fits')
  badpixfile=repstr(checkfile, 'FinalCR/'+snname+'_check.fits','badpix/'+snname+'.pix')
  print, checkfile
  print, fluxfile
  print, badpixfile
  ck=mrdfits(checkfile, 1)
  hand_mask=ck.indiv*!values.f_nan
  n_indiv = (size(ck.indiv, /dim))[0]
  infiles=strarr(n_indiv)
  for i=0, n_indiv-1 do begin
;     soplot, ck.indiv[i,*], ps=10, linestyle=1
     s='  '
     if i eq i_indiv then s='* '
     print, s+sxpar(headfits(fluxfile), 'INFILE'+strn(i,F='(I02)'))
     infiles[i] = sxpar(headfits(fluxfile), 'INFILE'+strn(i,F='(I02)'))
  endfor

  if (file_info(badpixfile)).exists then begin
     ;; open CR bad pixel file
     openr, lun, badpixfile, /get_lun
     ;; first line is the name of a single exposure file
     a=''
     readf, lun, a
     filename = a
     while ~eof(lun) do begin
        ;; next line of file
        a=''
        readf, lun, a
        if strpos(a, ' ') eq -1 then begin ;; found a filename
           ;; store previous exposure's data in structure
           badpixeldata0 = {filename:filename, badpixels:ptr_new(badpixels)}
           if n_elements(badpixeldata) eq 0 then begin
              badpixeldata = badpixeldata0
           endif else begin
              badpixeldata = [badpixeldata, badpixeldata0]
           endelse
           filename = a
           delvarx, badpixels
        endif else begin ;; found a pixel range
           range = fix(strsplit(a, ' ', /extract))
           if n_elements(badpixels) eq 0 then begin
              badpixels = indgen(range[1]-range[0]+1)+range[0]
           endif else begin
              badpixels = [badpixels, indgen(range[1]-range[0]+1)+range[0]]
           endelse
        endelse
     endwhile
     close, lun
     free_lun, lun
     ;; end of file encountered,
     ;; store data for last exposure
     badpixeldata0 = {filename:filename, badpixels:ptr_new(badpixels)}
     if n_elements(badpixeldata) eq 0 then begin
        badpixeldata = badpixeldata0
     endif else begin
        badpixeldata = [badpixeldata, badpixeldata0]
     endelse
     ;; now apply masked pixels to influx and inivar
     for i=0, n_elements(infiles)-1 do begin
        w=(where('Science_box/'+infiles[i] eq badpixeldata.filename))[0]
        if w eq -1 then continue
        hand_mask[i, *badpixeldata[w].badpixels] = 1
     endfor
  endif
  hand_mask[where(~finite(hand_mask))]=!values.f_nan
  splot, ck.coadd, ps=10
  for i=0, n_indiv-1 do begin
     soplot, ck.indiv[i,*], ps=10, linestyle=1
  endfor
  soplot, ck.indiv[i_indiv,*], ps=10, linestyle=2, color='magenta'
  soplot, ck.indiv[i_indiv,*]*hand_mask[i_indiv,*], ps=4, color='red'
  print, 'emacs '+badpixfile
end
