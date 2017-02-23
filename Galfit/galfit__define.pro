Function Galfit::Init, $
   params=params, $
   constraints=constraints, $
   components=components, $
   image=image, $
   errim=errim, $
   mskim=mskim, $
   psf=psf, $
   prefix=prefix, $
   outdir=outdir, $
   header=header
  if n_elements(params)      ne 0 then self.params = ptr_new(params)
  if n_elements(constraints) ne 0 then self.constraints = ptr_new(constraints)
  if n_elements(components)  ne 0 then self.components = ptr_new(components)
  if n_elements(image)       ne 0 then self.image = ptr_new(image)
  if n_elements(errim)       ne 0 then self.errim = ptr_new(errim)
  if n_elements(mskim)       ne 0 then self.mskim = ptr_new(mskim)
  if n_elements(psf)         gt 1 then self.psf = ptr_new(psf)
  if n_elements(prefix)      ne 0 then self.prefix = prefix
  if n_elements(outdir)      ne 0 then self.outdir = outdir $
  else begin
     cd, current=pwd
     self.outdir = pwd
  endelse
  self.outdir = directoryify(self.outdir)

  if n_elements(header)  ne 0 then self.header = ptr_new(header)

  return, 1
end

Pro Galfit::Execute, gzip=gzip
  self->Prepare
  file_delete, self.outdir+'run1/', /recursive, /allow_nonexistent
  spawn, 'cp -R '+self.outdir+'tmp/ '+self.outdir+'run1/'
  self->Run
  self->ReadFitLog
  if (file_info(self.outdir+'tmp/fit.log')).exists then begin  ;;if successful
;     spawn, 'mv '+self.outdir+'tmp/fit.log '+self.outdir+self.prefix+'.log'
;     spawn, 'mv '+self.outdir+'tmp/' $
;            +self.prefix+'.control ' $
;            +self.outdir+self.prefix+'.control'
;     spawn, 'mv '+self.outdir+'tmp/' $
;            +self.prefix+'.constraints ' $
;            +self.outdir+self.prefix+'.constraints'

     self.model = ptr_new(mrdfits( self.outdir+'tmp/imgblock.fits', 2, /silent ))
     self.resid = ptr_new(mrdfits( self.outdir+'tmp/imgblock.fits', 3, /silent ))

     self->RePrepare
     file_delete, self.outdir+'run2/', /recursive, /allow_nonexistent
     spawn, 'cp -R '+self.outdir+'tmp/ '+self.outdir+'run2/'
     self->Run
     self.intsb = ptr_new(mrdfits( self.outdir+'tmp/imgblock.fits', 3, /silent ))

     if ~ptr_valid(self.errim) then errim = *self.image*!values.f_nan $
     else errim = *self.errim
     struct = { image: *self.image, $
                errim: errim, $
                mskim: *self.mskim, $
                model: *self.model, $
                resid: *self.resid, $
                intsb: *self.intsb }

     results = (*self.components)[0]->Results() ;;get results for main galaxy
     struct = struct_addtags( struct, results )

     ;get sky if supplied
     comptypes = strarr(n_elements(*self.components))
     for i=0, n_elements(*self.components)-1 do begin
        comptypes[i] = obj_class((*self.components)[i])
     endfor
     wsky = where(comptypes eq 'GALFITSKY')
     if wsky[0] ne -1 then begin
        skycomponent = (*self.components)[wsky[0]]
        skyresults = skycomponent->Results()
        struct = struct_addtags( struct, skyresults )
     endif

     mwrfits, struct, self.outdir+self.prefix+'.fits', /create
     if keyword_set(gzip) then begin
        spawn, 'gzip -f '+self.outdir+self.prefix+'.fits'
     endif
     cd, current=pwd
     pwd += '/'
     cd, self.outdir
     spawn, 'tar -cvjpf '+self.prefix+'.bz2 run1 run2 > /dev/null'
     cd, pwd
  endif else self.fail = 1
  file_delete, self.outdir+'tmp/', /recursive, /allow_nonexistent
  file_delete, self.outdir+'run1/', /recursive, /allow_nonexistent
  file_delete, self.outdir+'run2/', /recursive, /allow_nonexistent
end

Pro Galfit::Run ;; actually step out of IDL to run galfit
  cd, current=pwd
  pwd += '/'
  cd, self.outdir+'tmp/'
  spawn, 'galfit '+self.prefix+'.control > /dev/null'
  cd, pwd
end

Pro Galfit::Prepare ;; prepare data/txt files for galfit
  file_mkdir, self.outdir
  file_mkdir, self.outdir+'tmp/'

  self->CreateDataFiles
  self->CreateConstraintFile
  self->CreateControlFile
end

Pro Galfit::RePrepare ;; change control file to not subtract target galaxy
  (*self.params).options = 2
  (*self.components)[0]->SetProperty, out=1
  self->CreateConstraintFile
  self->CreateControlFile
end

Function Galfit::Model
  file_mkdir, self.outdir
  file_mkdir, self.outdir+'tmp/'

  (*self.params).options = 2
  self->CreateDataFiles
  self->CreateControlFile
  self->CreateConstraintFile
  self->Run
  model = mrdfits( self.outdir+'tmp/imgblock.fits', 1, /silent ) $
          - mrdfits( self.outdir+'tmp/imgblock.fits', 2, /silent )
  file_delete, self.outdir+'tmp/', /recursive
  return, model
end

Function Galfit::Results
  if ~self.fail then begin
     results=(*self.components)[0]->Results()
     results = struct_addtags( results, {chi2:self.chi2, dof:self.dof, chi2perdof:self.chi2perdof} )
     return, results
  endif else return, -1
end

Pro Galfit::CreateDataFiles
  if ptr_valid(self.header) then $
     mwrfits, *self.image, self.outdir+'tmp/gal.fits', *self.header, /create $
  else $
     mwrfits, *self.image, self.outdir+'tmp/gal.fits', /create
  if ptr_valid(self.errim) then $
     mwrfits, *self.errim, self.outdir+'tmp/errim.fits', /create
  if ptr_valid(self.mskim) then $
     mwrfits, *self.mskim, self.outdir+'tmp/mskim.fits', /create
  if ptr_valid(self.psf) then $
     mwrfits, *self.psf, self.outdir+'tmp/psf.fits', /create
end

Pro Galfit::CreateConstraintFile
  if ~ptr_valid(self.constraints) then return
  openw, lun, self.outdir+'tmp/'+self.prefix+'.constraints', /get_lun
  for i=0, n_elements( *self.constraints ) - 1 do begin
     (*self.constraints)[i]->WriteToConstraintFile, lun
  endfor
  close, lun
  free_lun, lun
end

Pro Galfit::CreateControlFile
  openw, lun, self.outdir+'tmp/'+self.prefix+'.control', /get_lun
  self->WriteControlFileHeader, lun

  for i=0, n_elements( *self.components ) - 1 do begin
     (*self.components)[i]->WriteToControlFile, lun, i+1
  endfor

  close, lun
  free_lun, lun
end

Pro Galfit::WriteControlFileHeader, $
   lun

  printf, lun, '# IMAGE PARAMETERS'
  printf, lun, 'A) gal.fits            # Input data image (FITS file)'
  printf, lun, 'B) imgblock.fits       # Output data image block'
  if ~ptr_valid(self.errim) then $
     printf, lun, 'C) none                # Sigma image name (made from data if blank of "none")' $
  else $
     printf, lun, 'C) errim.fits          # Sigma image name (made from data if blank of "none")'
  if ptr_valid(self.psf) then begin
     printf, lun, 'D) psf.fits            # Input PSF image and (optional) diffusion kernel'
     printf, lun, 'E) '+string((*self.params).psffactor)+'                  # PSF fine sampling factor relative to data'
  endif else begin
     printf, lun, 'D)                     # Input PSF image and (optional) diffusion kernel'
     printf, lun, 'E)                     # PSF fine sampling factor relative to data'
  endelse

  if ~ptr_valid(self.mskim) then $
     printf, lun, 'F) none                # Bad pixel mask (FITS image or ASCII coord list)' $
  else $
     printf, lun, 'F) mskim.fits          # Bad pixel mask (FITS image or ASCII coord list)'

  if ~ptr_valid(self.constraints) then $
     printf, lun, 'G) none                # File with parameter constraints (ASCII file)' $
  else $
     printf, lun, 'G) '+self.prefix+'.constraints  # File with parameter constraints (ASCII file)'

  printf, lun, 'H) '+string((*self.params).fit_xmin)+string((*self.params).fit_xmax)+ $
          string((*self.params).fit_ymin)+string((*self.params).fit_ymax)+' #Image region to fit  (xmin xmax ymin ymax)'
  printf, lun, 'I) '+string((*self.params).conv_x)+string((*self.params).conv_y)+'  # Size of the convolution box (x y)'
  printf, lun, 'J) '+string((*self.params).zeropoint)+'       # Magnitude photometric zeropoint'
  printf, lun, 'K) '+string((*self.params).platescalex)+string((*self.params).platescaley)+'  # Plate scale (dx dy)  [arcsec per pixel]'
  printf, lun, 'O) regular             # Display type (regular, curses, both)'
  printf, lun, 'P) '+string((*self.params).options)+'        # Options: 0=normal run; 1,2=make model/imgblock & quit'
  printf, lun, 'S) 0                   # Modify/create objects interactively?'
  printf, lun
end

Pro Galfit::ReadFitLog
  if ~(file_info(self.outdir+'tmp/fit.log')).exists then return
  openr, lun, self.outdir+'tmp/fit.log', /get_lun
  a=''
  for i=0,6 do readf, lun, a
  objnum = 0
  while ~eof(lun) do begin
     readf, lun, a
     result = strsplit( a, ' :()[],=*', /extract )
     case result[0] of
        'sersic': begin
           (*self.components)[objnum]->SetProperty, xpos=double( result[1] ),  $
                                                    ypos=double( result[2] ), $
                                                    mag=double( result[3] ), $
                                                    re=double( result[4] ), $
                                                    n=double( result[5] ), $
                                                    boa=double( result[6] ), $
                                                    pa=double( result[7] )
           a=''
           readf, lun, a
           errors = strsplit( a, ' :()[],=*', /extract )
           (*self.components)[objnum]->SetProperty, errxpos=double( errors[0] ), $
                                                    errypos=double( errors[1] ), $
                                                    errmag=double( errors[2] ), $
                                                    errre=double( errors[3] ), $
                                                    errn=double( errors[4] ), $
                                                    errboa=double( errors[5] ), $
                                                    errpa=double( errors[6] )
           objnum += 1
        end
        'sky': begin
           (*self.components)[objnum]->SetProperty, sky=double( result[3] ), $
                                                    dskydx=double( result[4] ), $
                                                    dskydy=double( result[5] )
           a=''
           readf, lun, a
           errors = strsplit( a, ' :()[],=*', /extract )
           (*self.components)[objnum]->SetProperty, errsky=double( errors[0] ), $
                                                    errdskydx=double( errors[1] ), $
                                                    errdskydy=double( errors[2] )
           objnum += 1
        end
        'psf': begin
           (*self.components)[objnum]->SetProperty, xpos=double( result[1] ), $
                                                    ypos=double( result[2] ), $
                                                    mag=double( result[3] )
           a=''
           readf, lun, a
           errors = strsplit( a, ' :()[],=*', /extract )
           (*self.components)[objnum]->SetProperty, errxpos=double( errors[0] ), $
                                                    errypos=double( errors[1] ), $
                                                    errmag=double( errors[2] )
           objnum += 1
        end
        'Chi^2': begin
           self.chi2 = double( result[1] )
           self.dof = double( result[3] )
           a=''
           readf, lun, a
           result = strsplit( a, ' :()[],=*', /extract )
           self.chi2perdof = double( result[1] )
        end
        else :
     endcase
  endwhile
  close, lun
  free_lun, lun
end

Pro Galfit::CleanUp
  ptr_free, self.params
  ptr_free, self.image
  ptr_free, self.errim
  ptr_free, self.mskim
  ptr_free, self.model
  ptr_free, self.resid
  ptr_free, self.intsb
  if ptr_valid(self.psf) then ptr_free, self.psf
  if ptr_valid(self.components) then obj_destroy, *self.components
  if ptr_valid(self.constraints) then obj_destroy, *self.constraints
  ptr_free, self.constraints
  ptr_free, self.components
end

Pro Galfit__Define, struct
  struct = { Galfit, $
             Inherits JEMobject, $
             params:      ptr_new(), $
             constraints: ptr_new(), $
             components:  ptr_new(), $
             image:       ptr_new(), $
             errim:       ptr_new(), $
             mskim:       ptr_new(), $
             model:       ptr_new(), $
             resid:       ptr_new(), $
             intsb:       ptr_new(), $
             psf:         ptr_new(), $
             header:      ptr_new(), $
             outdir:      '', $
             prefix:      '', $
             chi2:        0.d, $
             dof:         0.d, $
             chi2perdof:  0.d, $
             fail:        0 }
end

