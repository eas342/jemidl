Function AddSpectra1, $
   dir, $      ;; directory in which spectra to add reside
   filename, $ ;; catalog name (in directory dir)
   sexcat      ;; sextractor catalog

  thisexcept = !except
  !except=0

  dir = directoryify(dir)
  SpecCatFile = dir+filename
  openr, lun, SpecCatFile, /get_lun
  f=''
  readf, lun, f, format='(A)' ;;get the source
  source=f
  readf, lun, mjd, format='(D)';;get the MJD
  a=''
  readf, lun, a, format='(A)' ;;get the x/y offsets. (in seconds of arc)
  a=strsplit( a, ' ', /extract )
  a=double(a)
  dra=a[0]
  ddec = a[1]
  a=''
  readf, lun, a, format='(A)' ;;get the first line
  a = strsplit( a, ' ', /extract )
  name    = a[0]
  rastr   = a[1]
  decstr  = a[2]
  z       = a[3]
  zerr    = a[4]
  zqual   = a[5]
  comment = a[6]
  while not eof(lun) do begin ;;get the rest of the lines
     a = ''
     readf, lun, a, format='(A)'
     a = strsplit( a,' ', /extract )
     if strmid(a[0],0,1) eq '#' then continue ;; skip comment lines
     name    = [name, a[0]]
     rastr   = [rastr, a[1]]
     decstr  = [decstr, a[2]]
     z       = [z, a[3]]
     zerr    = [zerr, a[4]]
     zqual   = [zqual, a[5]]
     comment = [comment, a[6]]
  endwhile
  close, lun
  free_lun, lun

  ngals = n_elements(name)
  ra = dblarr(ngals)
  dec = dblarr(ngals)
  shiftra = dblarr(ngals)
  shiftdec = dblarr(ngals)

  ;; get ra/dec in decimal degrees
  origra  = dblarr(ngals)
  origdec = dblarr(ngals)
  for i=0, ngals-1 do begin
     get_coords, coords, instring=rastr[i]+' '+decstr[i]
     ra[i]  = coords[0]*15.
     dec[i] = coords[1]
     origra[i] = ra[i]
     origdec[i] = dec[i]
     ra[i]  += dra/cos(dec[i]*!dpi/180.d)/3600. ;;arcsec -> degrees
     dec[i] += ddec/3600.
  endfor

  ;; now calculate shift
  for i=0, ngals-1 do begin
     dist2 = (sexcat.alphawin_j2000 - ra[i])^2*cos( dec[i]*!pi/180. )^2 $
             + (sexcat.deltawin_j2000 - dec[i])^2
     md2 = min( dist2, m )*3600.^2
     shiftra[i]  = ( ra[i]  - sexcat[m].alphawin_j2000 )*3600.
     shiftdec[i] = ( dec[i] - sexcat[m].deltawin_j2000 )*3600.
  endfor
  medshiftra  = median( shiftra[ where( abs( shiftra  ) lt 2. )] )/3600.d
  medshiftdec = median( shiftdec[where( abs( shiftdec ) lt 2. )] )/3600.d

  print
  print, 'medshiftra:', medshiftra*3600.
  print, 'medshiftdec:', medshiftdec*3600.
  print

  case source of
     'DEIMOS' : begin
        for i=0, ngals-1 do begin
           counter, i+1, ngals, 'Adding DEIMOS spectra: '
           if z[i] eq '?' or z[i] eq '??' then begin ;;targeted but failed to get redshift.
              z[i]       = -1.d
              zerr[i]    = -1.d
              zqual[i]   = 'F'
              comment[i] = ''
           endif

           dist2 = (sexcat.alphawin_j2000 - (ra[i] - medshiftra))^2*cos( dec[i]*!pi/180. )^2 $
                   + (sexcat.deltawin_j2000 - (dec[i] - medshiftdec))^2
           md2 = min( dist2, m )*3600.d^2

           deimos_getspec, dir+name[i], wave, flux, ivar
           if double(z[i]) eq -1.d then begin
              orig_z = obj_new( 'DopplerShift', z=0.d )
              new_z  = obj_new( 'DopplerShift', z=0.d )
           endif else begin
              orig_z = obj_new( 'DopplerShift', z=double(z[i]) )
              new_z  = obj_new( 'DopplerShift', z=double(z[i]) )
           endelse
           spectra1 = obj_new( 'spectrum', wave, flux, ivar=ivar, $
                               unit='flambda', orig_z=orig_z, new_z=new_z)
           if md2 lt 1.0^2 then begin ;; found in the ACS image
              clusterid = sexcat[m].clusterid
              galid = sexcat[m].galid
           endif else begin ;; not identified in the ACS image
              clusterid = sexcat[0].clusterid
              galid = -1
           endelse

           specstruct = { SpectrumTab, $
                          clusterid: clusterid, $
                          galid:     galid, $
                          source:    source, $
                          ra:        ra[i] - medshiftra, $
                          dec:       dec[i] - medshiftdec, $
                          MJD:       MJD, $
                          z:         double(z[i]), $
                          zerr:      double(zerr[i]), $
                          zqual:     zqual[i], $
                          comment:   comment[i], $
                          o2_EW:     !values.f_nan, $
                          o2_EWerr:  !values.f_nan, $
                          filename:  dir+name[i], $
                          spectrum:  spectra1 }
           if n_elements(spectra) eq 0 then spectra=specstruct $
           else spectra=[spectra, specstruct]
        endfor
     end ;;'DEIMOS'
     'FORS1/2' : begin
        for i=0, n_elements(name)-1 do begin
           counter, i+1, ngals, 'Adding FORS1/2 spectra: '
           if z[i] eq '?' or z[i] eq '??' then begin
              z[i]       = -1.d
              zerr[i]    = -1.d
              zqual[i]   = 'F'
              comment[i] = ''
           endif

           dist2 = (sexcat.alphawin_j2000 - ra[i] + medshiftra)^2*cos( dec[i]*!pi/180. )^2 $
                   + (sexcat.deltawin_j2000 - dec[i] + medshiftdec)^2
           md2 = min( dist2, m )*3600.d^2
           
           files=file_search( dir+'**/*'+name[i]+'*' )
           werr=where( strpos( files, 'err' ) ne -1, complement=wspec )
           if werr[0] eq -1 then continue
           specfile= files[wspec]
           errfile = files[werr]
           flux = mrdfits( specfile[0], /silent )
           err  = mrdfits( errfile[0],  /silent )
           if n_elements(flux) ne n_elements(err) then continue
           h=headfits(specfile[0])
           cd1_1  = double( sxpar( h, 'CD1_1'  ) )
           crval1 = double( sxpar( h, 'CRVAL1' ) )
           if crval1 eq 0 then begin
              string = sxpar( h, 'WAT2_001' )
              stuff = strsplit( string, ' ', /extract)
              crval1 = double( stuff[6] ) 
              cd1_1 = double( stuff[7] ) 
           endif
           wave = findgen( n_elements(flux) )*cd1_1+crval1
           if double(z[i]) eq -1.d then begin
              orig_z = obj_new( 'DopplerShift', z=0.d )
              new_z  = obj_new( 'DopplerShift', z=0.d )
           endif else begin
              orig_z = obj_new( 'DopplerShift', z=double(z[i]) )
              new_z  = obj_new( 'DopplerShift', z=double(z[i]) )
           endelse
           spectra1 = obj_new( 'spectrum', wave, flux, ivar=1./err^2, $
                               unit='flambda', orig_z=orig_z, new_z=new_z)
           if md2 lt 1.0^2 then begin ;; found in the ACS image
              clusterid = sexcat[m].clusterid
              galid = sexcat[m].galid
           endif else begin ;; not identified in the ACS image
              clusterid = sexcat[0].clusterid
              galid = -1
           endelse

           specstruct = { SpectrumTab, $
                          clusterid: clusterid, $
                          galid:     galid, $
                          source:    source, $
                          ra:        ra[i] - medshiftra, $
                          dec:       dec[i] - medshiftdec, $
                          MJD:       MJD, $
                          z:         double(z[i]), $
                          zerr:      double(zerr[i]), $
                          zqual:     zqual[i], $
                          comment:   comment[i], $
                          o2_EW:     !values.f_nan, $
                          o2_EWerr:  !values.f_nan, $
                          filename:  dir+name[i], $
                          spectrum:  spectra1 }
           if n_elements(spectra) eq 0 then spectra=specstruct $
           else spectra=[spectra, specstruct]
        endfor
     end ;;'FORS1/2'
     'FOCAS' : begin
        for i=0, n_elements(name)-1 do begin
           counter, i+1, ngals, 'Adding FOCAS spectra: '
           if z[i] eq '?' then begin
              z[i]       = -1.d
              zerr[i]    = -1.d
              zqual[i]   = 'F'
              comment[i] = ''
           endif

           dist2 = (sexcat.alphawin_j2000 - ra[i] + medshiftra)^2*cos( dec[i]*!pi/180. )^2 $
                   + (sexcat.deltawin_j2000 - dec[i] + medshiftdec)^2
           md2 = min( dist2, m )*3600.d^2
           
           readcol, dir+name[i]+'.txt', wave, flux, format='D,D', /silent
           readcol, dir+name[i]+'error.txt', wave, err, format='D,D', /silent
           ivar = 1.d/err^2
           w = where( 1-finite(ivar) )
           if w[0] ne -1 then ivar[w]=0.d
           
           if double(z[i]) eq -1.d then begin
              orig_z = obj_new( 'DopplerShift', z=0.d )
              new_z  = obj_new( 'DopplerShift', z=0.d )
           endif else begin
              orig_z = obj_new( 'DopplerShift', z=double(z[i]) )
              new_z  = obj_new( 'DopplerShift', z=double(z[i]) )
           endelse
           spectra1 = obj_new( 'spectrum', wave, flux, ivar=ivar, $
                               unit='flambda', orig_z=orig_z, new_z=new_z )
           if md2 lt 1.0^2 then begin ;; found in the ACS image
              clusterid = sexcat[m].clusterid
              galid = sexcat[m].galid
           endif else begin ;; not identified in the ACS image
              clusterid = sexcat[0].clusterid
              galid = -1
           endelse

           specstruct = { SpectrumTab, $
                          clusterid: clusterid, $
                          galid:     galid, $
                          source:    source, $
                          ra:        ra[i] - medshiftra, $
                          dec:       dec[i] - medshiftdec, $
                          MJD:       MJD, $
                          z:         double(z[i]), $
                          zerr:      double(zerr[i]), $
                          zqual:     zqual[i], $
                          comment:   comment[i], $
                          o2_EW:     !values.f_nan, $
                          o2_EWerr:  !values.f_nan, $
                          filename:  dir+name[i], $
                          spectrum:  spectra1 }
           if n_elements(spectra) eq 0 then spectra=specstruct $
           else spectra=[spectra, specstruct]
        endfor
     end ;;'FOCAS'
     else : begin ;; 'Literature sources'
        for i=0, n_elements(name)-1 do begin
           counter, i+1, ngals, 'Adding Literature spectra: '           
           dist2 = (sexcat.alphawin_j2000 - ra[i] + medshiftra)^2*cos( dec[i]*!pi/180. )^2 $
                   + (sexcat.deltawin_j2000 - dec[i] + medshiftdec)^2
           md2 = min( dist2, m )*3600.d^2

           if md2 lt 1.0^2 then begin ;;in the ACS image
              clusterid = sexcat[m].clusterid
              galid = sexcat[m].galid
           endif else begin
              clusterid = sexcat[0].clusterid
              galid = -1
           endelse
           specstruct = { SpectrumTab, $
                          clusterid: clusterid, $
                          galid:     galid, $
                          source:    source, $
                          ra:        ra[i] - medshiftra, $
                          dec:       dec[i] - medshiftdec, $
                          MJD:       MJD, $
                          z:         double(z[i]), $
                          zerr:      double(zerr[i]), $
                          zqual:     zqual[i], $
                          comment:   comment[i], $
                          o2_EW:     !values.f_nan, $
                          o2_EWerr:  !values.f_nan, $
                          filename:  dir+name[i], $
                          spectrum:  obj_new() }
           if n_elements(spectra) eq 0 then spectra=specstruct $
           else spectra=[spectra, specstruct]
        endfor
     end ;;'Literature'
  endcase
  blah = check_math()
  !except=thisexcept
  return, spectra
end
