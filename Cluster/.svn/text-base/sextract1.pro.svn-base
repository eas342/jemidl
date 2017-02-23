Pro sextract1, $
   sexfile, $
   paramfile, $
   imagefile, $
   zeropoint, $
   outname, $
   weightfile=weightfile, $
   weighttype=weighttype, $
   segment=segment, $
   background=background

  nfiles=n_elements(imagefile)
  if nfiles eq 1 then begin ;;single image mode
     hdr = headfits(imagefile)
     exptime = sxpar( hdr, 'EXPTIME' )
     cmd = 'sex '+imagefile $
           + ' -c '+sexfile $
           + ' -PARAMETERS_NAME '+paramfile $
           + ' -MAG_ZEROPOINT '+strtrim(zeropoint,2) $
           + ' -GAIN '+str(exptime) $
           + ' -CATALOG_NAME '+outname
     if n_elements(weightfile) ne 0 then begin ;;add optional weight image
        cmd += ' -WEIGHT_IMAGE '+weightfile $
               + ' -WEIGHT_TYPE '+weighttype
     endif
     if n_elements(segment) ne 0 then begin ;;optional output: segmentation image
        checkimage_type = 'SEGMENTATION'
        checkimage_name = segment
     endif
     if n_elements(background) ne 0 then begin ;;optional output: background image
        if n_elements(checkimage_type) eq 0 then checkimage_type = 'BACKGROUND' $
        else checkimage_type += ',BACKGROUND'
        if n_elements(checkimage_name) eq 0 then checkimage_name = background $
        else checkimage_name += ','+str(background)
     endif
     if n_elements(checkimage_type) ne 0 then begin
        cmd += ' -CHECKIMAGE_TYPE '+checkimage_type $
               + ' -CHECKIMAGE_NAME '+checkimage_name
     endif
  endif else begin ;;dual image mode
     hdr = headfits(imagefile[1])
     exptime = sxpar( hdr, 'EXPTIME' )
     cmd = 'sex '+imagefile[0]+','+imagefile[1] $
           + ' -c '+sexfile $
           + ' -PARAMETERS_NAME '+paramfile $
           + ' -MAG_ZEROPOINT '+str(zeropoint[1]) $
           + ' -GAIN '+str(exptime) $
           + ' -CATALOG_NAME ' + outname
     if n_elements(weightfile[0]) ne 0 then begin ;;optional input: weight image
        cmd += ' -WEIGHT_IMAGE '+weightfile[0]+','+weightfile[1] $
               + ' -WEIGHT_TYPE '+weighttype
     endif
  endelse
  print, cmd
  spawn, cmd
  if n_elements(segment) ne 0 and nfiles eq 1 then begin
     spawn, 'gzip -f '+segment
  endif
end
