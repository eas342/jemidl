;;;;;;;;;;;;;;;;;;;;;;
;;  Init
;;Initialize object
Function Cluster3::Init, $
   _extra=extra

  self->SetProperty, _extra=extra
  return, 1
end


;;;;;;;;;;;;;;;;;;;;;;
;;  PrepareSExtractorFiles
;;Prepare SExtractor input files
Pro Cluster3::PrepareSExtractorFiles, band

  if n_elements(band) eq 0 then band = (*self.band)
  scratchdir = '/scratch-local/jmeyers314/sex'+self.clusterid+'/'
  file_mkdir, scratchdir
  
  for iband=0, n_elements(band)-1 do begin
     xband = (where(band[iband] eq *self.band))[0]
     split = strsplit( (*self.filehead)[xband], '/', /extract )
     name = split[n_elements(split)-1]
     cmd = 'gunzip -c ' $
           +self.srcdir+(*self.filehead)[xband] $
           +'_SCI.fits.gz > ' $
           +scratchdir+name $
           +'_SCI.fits'
     spawn, cmd
     cmd = 'gunzip -c ' $
           +self.srcdir+(*self.filehead)[xband] $
           +'_ERR.fits.gz > ' $
           +scratchdir+name $
           +'_ERR.fits'
     spawn, cmd
  endfor
end


;;;;;;;;;;;;;;;;;;;;;;
;;  CleanUpSExtractorFiles
;;Cleanup SExtractor input files
Pro Cluster3::CleanUpSExtractorFiles
  scratchdir='/scratch-local/jmeyers314/sex'+self.clusterid+'/'
  cmd = 'rm -rf '+scratchdir
  spawn, cmd
end


;;;;;;;;;;;;;;;;;;;;;;
;;  SExtract
;;Run SExtractor for (optionally) multiple bands.
Function Cluster3::SExtract, SExfile, $
                             paramfile, $
                             band, $
                             outdir=outdir, $
                             noprep=noprep, $
                             _extra=extra
;; _extra: segment=segment,
;;         background=background
;;         these last two only effectively apply to the index-image


  if n_elements(outdir) eq 0 then outdir=self.rootdir
  outdir=directoryify(outdir)
  file_mkdir, outdir
  ;;determine necessary file preparation
  if n_elements(band) eq 0 then begin
     band = (*self.band)
     prepband = (*self.band)
  endif else begin
     prepband = band
     if ~member( (*self.band)[0], prepband ) $
     then prepband=[(*self.band)[0], prepband]
  endelse
  if ~keyword_set(noprep) then self->PrepareSExtractorFiles, prepband
  scratchdir='/scratch-local/jmeyers314/sex'+self.clusterid+'/'

  for iband=0, n_elements(band)-1 do begin
     xband = (where(band[iband] eq *self.band))[0]
     split = strsplit( (*self.filehead)[xband], '/', /extract )
     name = split[n_elements(split)-1]
     split0 = strsplit( (*self.filehead)[0], '/', /extract )
     name0 = split0[n_elements(split0)-1]

     ;;SExtractor argument construction
     if xband eq 0 then begin ;;index image: use single image mode
        imagefile = scratchdir+name+'_SCI.fits'
        weightfile = scratchdir+name+'_ERR.fits'
        zeropoint = (*self.zeropoint)[0]
     endif else begin ;;non-index image: use dual image mode
        imagefile = scratchdir+[name0, name]+'_SCI.fits'
        weightfile = scratchdir+[name0, name]+'_ERR.fits'
        zeropoint = [(*self.zeropoint)[0], (*self.zeropoint)[xband]]
     endelse

     ;;run SExtractor
     SExtract1, SExfile, paramfile, $
                imagefile, zeropoint, outdir+band[iband]+'cat.fits', $
                weightfile=weightfile, weighttype='MAP_RMS', _extra=extra
     
     ;;import results
     cat0 = mrdfits( outdir+band[iband]+'cat.fits', 2 )
     
     ;;add some tags to the imported catalog
     nobjs = (size( cat0, /dim ))[0]
     s = { clusterid:self.clusterid, galid:0L, band:band[iband] }
     s = replicate( s, nobjs )
     s.galid=lindgen(nobjs)
     cat0 = struct_addtags( s, struct_trimtags( cat0, except_tags='NUMBER' ) )
     if n_elements(cat) eq 0 then cat=cat0 $
     else cat=[[cat], [cat0]]
  endfor

  ;;cleanup files
  if ~keyword_set(noprep) then self->CleanUpSExtractorFiles
  return, cat
end


;;;;;;;;;;;;;;;;;;;;;;
;;  SExtractColdHot
;;Do a double SExtraction for (optionally) multiple bands
Function Cluster3::SExtractColdHot, coldSExfile, $
                                    hotSExfile, $
                                    paramfile, $
                                    band, $
                                    outdir=outdir, $
                                    noprep=noprep, $
                                    store=store

  if n_elements(outdir) eq 0 then outdir=self.rootdir
  outdir=directoryify(outdir)

  ;;determine necessary file preparation
  if n_elements(band) eq 0 then begin
     band = (*self.band)
     prepband = band
  endif else begin
     prepband = band
     if ~member( (*self.band)[0], prepband ) $
     then prepband=[(*self.band)[0], prepband]
  endelse
  if ~keyword_set(noprep) then self->PrepareSExtractorFiles, prepband
  nbands = n_elements(band)

  ;;start with the cold SExtraction
  coldcat = self->SExtract( coldSExfile, paramfile, band, /noprep, $
                            outdir=outdir, segment=outdir+'cold_segment.fits' )
  ncoldobjs = (size(coldcat, /dim))[0]
  coldcat = struct_addtags( coldcat, replicate({cold:1}, ncoldobjs, nbands) )
  for iband=0, n_elements(band)-1 do begin
     filehead = outdir+(*self.band)[iband]+'cat'
     file_move, filehead+'.fits', filehead+'_cold.fits', /overwrite
  endfor

  ;;do the hot SExtraction
  hotcat = self->SExtract( hotSExfile, paramfile, band, /noprep, $
                           outdir=outdir, segment=outdir+'hot_segment.fits' )
  nhotobjs = (size( hotcat, /dim ))[0]
  hotcat = struct_addtags( hotcat, replicate({cold:0}, nhotobjs, nbands) )
  for iband=0, n_elements(band)-1 do begin
     filehead = outdir+(*self.band)[iband]+'cat'
     file_move, filehead+'.fits', filehead+'_hot.fits', /overwrite
  endfor
  if ~keyword_set(noprep) then self->PrepareSExtractorFiles, prepband

  ;;read in cold segmentation map and expand it 2 pixels
  segment = mrdfits(outdir+'cold_segment.fits.gz')
  sz=size( segment, /dim )
  xsize=sz[0]
  ysize=sz[1]
  segment = segment ne 0
  for ix=-2,2 do begin
     segment = segment or shift(segment, ix, 0)
  endfor
  for iy=-2,2 do begin
     segment = segment or shift(segment, 0, iy)
  endfor

  ;;now disregard hot objects already detected in the cold SExtraction
  cat = coldcat
  for ihot=0, nhotobjs-1 do begin
     counter, ihot+1, nhotobjs, "merging catalogs: "
     hotcat1 = hotcat[ihot,*]
     x=floor(hotcat1[0].xwin_image)-1
     y=floor(hotcat1[0].ywin_image)-1
     if x lt 0 or x ge xsize or y lt 0 or y ge ysize then begin
        cat = [cat, hotcat1]
        continue
     endif
     if segment[x,y] eq 0 then cat = [cat, hotcat1]
  endfor
  print
  nobjs = (size( cat, /dim ))[0]
  cat.galid = rebin( lindgen(nobjs), nobjs, nbands )
  if keyword_set(store) then begin
     ptr_free, self.SExcat
     self.SExcat = ptr_new(cat)
  endif
  if ~keyword_set(noprep) then self->CleanUpSExtractorFiles
  return, cat
end


;;;;;;;;;;;;;;;;;;;;;;
;;  Galapagos
;;Algorithm to setup massively parallel galfiting
Function Cluster3::Galapagos, store=store
  if ~ptr_valid(self.SExcat) then begin
     message, 'Cannot run galapagos without first storing a SExtractor catalog', $
              /continue
     return, 0
  endif
  galapagoscat = galapagos1( *self.SExcat )
  if keyword_set(store) then begin
     ptr_free, self.galapagoscat
     self.galapagoscat = ptr_new(galapagoscat)
  endif
  return, galapagoscat
end


;;;;;;;;;;;;;;;;;;;;;;
;;  LoadImage
;;
Pro Cluster3::LoadImage, band=band, $
                         type=type
  nbands=n_elements(*self.band)
  if n_elements(band) eq 0 then band=(*self.band)
  if n_elements(type) eq 0 then type=['SCI','WHT','ERR','SKY','MSK','CVL','NXP']
  id = self.clusterid

  for itype=0, n_elements(type)-1 do begin
     if type[itype] eq 'SEGHOT' then begin
        self.seghot = ptr_new(mrdfits(self.rootdir+'hot_segment.fits.gz'))
        continue
     endif
     if type[itype] eq 'SEGCLD' then begin
        self.segcld = ptr_new(mrdfits(self.rootdir+'cold_segment.fits.gz'))
        continue
     endif
     for iband=0, n_elements(band)-1 do begin
        string1 = "if ~ptr_valid(self." $
                  +type[itype] $
                  +"im) then self." $
                  +type[itype] $
                  +"im = ptr_new(ptrarr(nbands))"
        junk = execute(string1)
        xband = (where(band[iband] eq *self.band))[0]
        string2 = "(*self." $
                  +type[itype] $
                  +"im)[xband] = ptr_new(mrdfits(self.srcdir+(*self.filehead)[xband]+'_" $
                  +type[itype] $
                  +".fits.gz'))"
        junk = execute(string2)
     endfor
  endfor
end


;;;;;;;;;;;;;;;;;;;;;;
;;  FreeImages
;;Cleanup image pointers
Pro Cluster3::FreeImages, type=type, band=band
  if n_elements(type) eq 0 $
  then type=['SCI','WHT','ERR','SKY','MSK','CVL','NXP','SEGHOT','SEGCLD']
  if n_elements(band) eq 0 then band=(*self.band)

  for itype=0, n_elements(type)-1 do begin
     if type[itype] eq 'SEGHOT' then begin
        ptr_free, self.seghot
        continue
     endif
     if type[itype] eq 'SEGCLD' then begin
        ptr_free, self.segcld
        continue
     endif
     string1 = "if ptr_valid(self." $
               +type[itype]+"im) then ptr_free, *self." $
               +type[itype]+"im"
     junk = execute(string1)
  endfor
end


;;;;;;;;;;;;;;;;;;;;;;
;;  LoadSNeFree
;;
Pro Cluster3::LoadSNeFree, snname
  sncoord = round(self->sncoord( snname, name1=name1 ))
  dir = '/home/scpdata03/joshimages/SNe/'
  num = fix(strmid( name1, 6 ))
  id = self.clusterid
  if id eq 'D' then id = 'O'
  subdir = 'CL-'+id+'-'+str(num, format='(I03)')+'/'
  filehead = (file_search(dir+subdir+'*_SCI.fits.gz'))[0]
  filehead = strmid( filehead, 0, strlen(filehead)-12 )
  snsci = mrdfits( filehead+'_SCI.fits.gz' , /silent )
  snerr = mrdfits( filehead+'_ERR.fits.gz' , /silent )
  snwht = mrdfits( filehead+'_WHT.fits.gz' , /silent )
  snnxp = mrdfits( filehead+'_NXP.fits.gz' , /silent )
  sci = *(*self.sciim)[0]
  err = *(*self.errim)[0]
  wht = *(*self.whtim)[0]
  nxp = *(*self.nxpim)[0]
  ptr_free, (*self.sciim)[0]
  ptr_free, (*self.errim)[0]
  ptr_free, (*self.whtim)[0]
  ptr_free, (*self.nxpim)[0]
  xmin=sncoord[0]-250
  xmax=sncoord[0]+250
  ymin=sncoord[1]-250
  ymax=sncoord[1]+250
  sci[xmin:xmax, ymin:ymax] = snsci
  err[xmin:xmax, ymin:ymax] = snerr
  wht[xmin:xmax, ymin:ymax] = snwht
  nxp[xmin:xmax, ymin:ymax] = snnxp
  (*self.sciim)[0] = ptr_new(sci)
  (*self.errim)[0] = ptr_new(err)
  (*self.whtim)[0] = ptr_new(wht)
  (*self.nxpim)[0] = ptr_new(nxp)
end

;;;;;;;;;;;;;;;;;;;;;;
;;  LoadGOODSSN
;;
Pro Cluster3::LoadGOODSSN, snname
  sncoord = round(self->sncoord(snname, name1=name1))
  dir = '/home/scpdata03/joshimages/SNe/'
  subdir = name1+'/'
  filehead = (file_search(dir+subdir+'*_SCI.fits.gz'))[0]
  filehead = strmid( filehead, 0, strlen(filehead)-12 )
  snsci = mrdfits( filehead+'_SCI.fits.gz' , /silent )
  snerr = mrdfits( filehead+'_ERR.fits.gz' , /silent )
  snwht = mrdfits( filehead+'_WHT.fits.gz' , /silent )
  snnxp = mrdfits( filehead+'_NXP.fits.gz' , /silent )
  if size( snsci, /n_dim ) eq 0 then return
  sci = *(*self.sciim)[0]
  err = *(*self.errim)[0]
  wht = *(*self.whtim)[0]
  nxp = *(*self.nxpim)[0]
  ptr_free, (*self.sciim)[0]
  ptr_free, (*self.errim)[0]
  ptr_free, (*self.whtim)[0]
  ptr_free, (*self.nxpim)[0]
  xmin=sncoord[0]-250
  xmax=sncoord[0]+250
  ymin=sncoord[1]-250
  ymax=sncoord[1]+250
  sci[xmin:xmax, ymin:ymax] = snsci
  err[xmin:xmax, ymin:ymax] = snerr
  wht[xmin:xmax, ymin:ymax] = snwht
  nxp[xmin:xmax, ymin:ymax] = snnxp
  (*self.sciim)[0] = ptr_new(sci)
  (*self.errim)[0] = ptr_new(err)
  (*self.whtim)[0] = ptr_new(wht)
  (*self.nxpim)[0] = ptr_new(nxp)
end

;;;;;;;;;;;;;;;;;;;;;;
;;  LoadHDR
;;
Pro Cluster3::LoadHDR
  self.hdr = ptr_new(headfits(self.srcdir+(*self.filehead)[0]+'_SCI.fits.gz'))
end

;;;;;;;;;;;;;;;;;;;;;;
;;  ImageSection
;;
Function Cluster3::ImageSection, p, $      ;; position to extract
                                 band=band, $
                                 type=type
  
  if n_elements(band) eq 0 then band=(*self.band)[0]
  if n_elements(type) eq 0 then type='SCI'
  xband = (where(band eq *self.band))[0]
  if type eq 'SEGHOT' then begin
     out=imagesection( *self.seghot, p )
  endif else if type eq 'SEGCLD' then begin
     out=imagesection( *self.segcld, p )
  endif else if type eq 'MSK' then begin
     string1 = "out=imagesection(*(*self." $
               +type+"im)[xband], p, background=1 )"
     junk = execute(string1)
  endif else begin
     string1 = "out=imagesection(*(*self." $
               +type+"im)[xband], p )"
     junk = execute(string1)
  endelse
  return, out
end


;;;;;;;;;;;;;;;;;;;;;;
;;  GalImage
;;
Function Cluster3::GalImage, gal, $
                             xsize=xsize, $
                             ysize=ysize, $
                             skysub=skysub, $
                             band=band, $
                             type=type
  if n_elements(band) eq 0 then band=(*self.band)[0]
  if n_elements(type) eq 0 then type='SCI'

  if n_elements(xsize) eq 0 or n_elements(ysize) eq 0 then begin
     gal1=(*self.galapagoscat)[gal]
     p=[gal1.xmin,gal1.ymin,gal1.xmax,gal1.ymax]
  endif else begin
     s=(*self.SExcat)[gal,0]
     x=s.xwin_image-1
     y=s.ywin_image-1
     p=floor([x-xsize/2, y-ysize/2, x+xsize/2, y+ysize/2])
  endelse
  img = self->ImageSection(p, band=band, type=type)
  if keyword_set(skysub) then begin
     skyvalue = (*self.skycat)[gal,where(band eq *self.band)].sky_value
     if finite(skyvalue) then img -= skyvalue
  endif
  return, img
end


;;;;;;;;;;;;;;;;;;;;;;
;;  Color
;;
Function Cluster3::Color, gals=gals, $
                          ap_radius=ap_radius, $
                          ap_factor=ap_factor, $
                          min_radius=min_radius, $
                          calc=calc, $
                          clean=clean, $
                          xpsf=xpsf, $
                          band1=band1, $
                          band2=band2, $
                          mags=mags, $
                          silent=silent, $
                          _extra=extra

  if n_elements(band1) eq 0 then band1 = 'i'
  if n_elements(band2) eq 0 then band2 = 'z'
  if n_elements(min_radius) eq 0 $
     and n_elements(ap_radius) eq 0 then begin
     min_radius=3.
     if ~keyword_set(silent) then begin
        print, 'WARNING: using minimum radius of 3.'
     endif
  endif
  if n_elements(ap_radius) ne 0 then min_radius=0.
  if n_elements(gals) eq 0 then gals=indgen( (size( *self.SExcat, /dim ))[0] )
  if n_elements(ap_factor) eq 0 and n_elements(ap_radius) eq 0 then ap_factor=1.
  clean = keyword_set(clean)
  xpsf = keyword_set(xpsf)

  iband1 = (where(band1 eq (*self.band)))[0]
  iband2 = (where(band2 eq (*self.band)))[0]
  color=fltarr( n_elements(gals) )
  if arg_present(mags) then mags = fltarr( n_elements(gals), 2 )
  diameters=[1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0, $
             11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0, $
             22.0,24.0,26.0,28.0,30.0,32.0,40.0,60.0,80.0,120.0,160.0,220.0]
  facs = [0.5, 0.75, 1.0, 1.25, 1.5]
  aps  = [3.0, 5.0, 10.0, 15.0, 20.0]

  if keyword_set(calc) then begin
     for igal=0, n_elements(gals)-1 do begin
        if n_elements(ap_factor) ne 0 then begin
           ;; correction for PSF convolved images
           if keyword_set(xpsf) then r = ap_factor*((*self.galfitcat)[gals[igal]].re+1.2) $
           else r = ap_factor*(*self.galfitcat)[gals[igal]].re
        endif else r = ap_radius
        if r lt 0 then begin
           color[igal]=!values.f_nan
           continue
        endif
        r = r > min_radius
        mag1 = self->ApPhot( r, band=band1, gals=gals[igal], $
                             clean=clean, xpsf=xpsf, _extra=extra, /silent )
        mag2 = self->ApPhot( r, band=band2, gals=gals[igal], $
                             clean=clean, xpsf=xpsf, _extra=extra, /silent )
        color[igal] = mag1 - mag2
        if arg_present(mags) then begin
           mags[igal, 0] = mag1
           mags[igal, 1] = mag2
        endif
     endfor
  endif else begin ;;don't explicity calculate
     if keyword_set(clean) then begin
        for igal=0, n_elements(gals)-1 do begin
           if n_elements(ap_factor) ne 0 then begin
              re = (*self.galfitcat)[gals[igal]].re
              mag1 = interpol( (*self.cleancat)[gals[igal],iband1].cleanfac, $
                               facs, ap_factor > (min_radius/re))
              mag2 = interpol( (*self.cleancat)[gals[igal],iband2].cleanfac, $
                               facs, ap_factor > (min_radius/re))
              color[igal] = mag1 - mag2
              if arg_present(mags) then begin
                 mags[igal, 0] = mag1
                 mags[igal, 1] = mag2
              endif
           endif else begin 
              mag1 = interpol( (*self.cleancat)[gals[igal],iband1].cleanap, $
                               aps, ap_radius > min_radius )
              mag2 = interpol( (*self.cleancat)[gals[igal],iband2].cleanap, $
                               aps, ap_radius > min_radius )
              color[igal] = mag1 - mag2
              if arg_present(mags) then begin
                 mags[igal, 0] = mag1
                 mags[igal, 1] = mag2
              endif
           endelse
        endfor
     endif else if keyword_set(xpsf) then begin
        for igal=0, n_elements(gals)-1 do begin
           if n_elements(ap_factor) ne 0 then begin
              re = (*self.galfitcat)[gals[igal]].re+1.2 ;;correction for PSF convolution
              mag1 = interpol( (*self.xpsfcat)[gals[igal],iband1].xpsffac, $
                               facs, ap_factor > (min_radius/re))
              mag2 = interpol( (*self.xpsfcat)[gals[igal],iband2].xpsffac, $
                               facs, ap_factor > (min_radius/re))
              color[igal] = mag1 - mag2
              if arg_present(mags) then begin
                 mags[igal, 0] = mag1
                 mags[igal, 1] = mag2
              endif
           endif else begin
              mag1 = interpol( (*self.xpsfcat)[gals[igal],iband1].xpsfap, $
                               aps, ap_radius > min_radius )
              mag2 = interpol( (*self.xpsfcat)[gals[igal],iband2].xpsfap, $
                               aps, ap_radius > min_radius )
              color[igal] = mag1 - mag2
              if arg_present(mags) then begin
                 mags[igal, 0] = mag1
                 mags[igal, 1] = mag2
              endif
           endelse
        endfor
     endif else begin
        for igal=0, n_elements(gals)-1 do begin ;;not PSF correction
           if n_elements(ap_factor) ne 0 then $
              radius = ap_factor*(*self.galfitcat)[gals[igal]].re $
           else radius = ap_radius
           radius = radius > min_radius
           mag1 = interpol((*self.SExcat)[gals[igal],iband1].mag_aper, $
                           diameters/2.0, radius)
           mag2 = interpol((*self.SExcat)[gals[igal],iband2].mag_aper, $
                           diameters/2.0, radius)
           color[igal] = mag1 - mag2
           if arg_present(mags) then begin
              mags[igal, 0] = mag1
              mags[igal, 1] = mag2
           endif
        endfor
     endelse
  endelse
  return, color
end

;;;;;;;;;;;;;;;;;;;;;;
;;  ApPhot
;;
Function Cluster3::ApPhot, radius, $
                           band=band, $
                           gals=gals, $
                           clean=clean, $
                           xpsf=xpsf, $
                           silent=silent, $
                           err=err, $
                           _extra=extra

  if n_elements(band) eq 0 then band=(*self.band)[0]
  iband = (where(band eq (*self.band)))[0]
  if n_elements(cband) eq 0 then cband=(*self.band)[1-iband]
  xband = (where(cband eq (*self.band)))[0]
  if n_elements(gals) eq 0 then gals=lindgen((size(*self.SExcat,/dim))[0])
  fluxes = fltarr(n_elements(gals),n_elements(radius))*!values.f_nan
  if arg_present(err) then $
     errflux = fltarr(n_elements(gals),n_elements(radius))*!values.f_nan

  if keyword_set(clean) then begin
     psf = *(*self.psf)[iband]
     psf = undersampleimage( psf, self.samplefactor )
  endif else if keyword_set(xpsf) then begin
     xpsf = undersampleimage( *(*self.psf)[xband], self.samplefactor )
     xpsf /= total(xpsf)
  endif

  for igal=0, n_elements(gals)-1 do begin
     if ~keyword_set(silent) then $
        counter, igal+1, n_elements(gals), 'Calculating aperture photometry: '
     SExcat1 = (*self.SExcat)[gals[igal],0]
     p = round([ SExcat1.xwin_image-1, $ ;;SEx coords to IDL coords
                 SExcat1.ywin_image-1, $
                 SExcat1.xwin_image-1, $
                 SExcat1.ywin_image-1 ]) + [-1,-1,1,1]*round((max(radius)*1.2) + 20)
     if total(p lt 0) ne 0 then begin
        fluxes[igal,*]=!values.d_nan
        continue
     endif
     subimage = self->ImageSection( p, band=band, type='SCI' )
     if arg_present(err) then err = self->ImageSection( p, band=band, type='ERR' )
     if ptr_valid(self.skycat) then sky = (*self.skycat)[gals[igal],iband].sky_value $
     else sky = 0.
     if ~finite(sky) then continue
     subimage -= sky
     if keyword_set(clean) then begin
        subimage = cleanimage( subimage, psf, 1.5 )
     endif else if keyword_set(xpsf) then begin
        subimage = convol( subimage, xpsf, /edge_truncate )
     endif
     fluxes[igal,*]=apphot( subimage, $
                            (SExcat1.xwin_image-0.5)-p[0], $  ;;SEx to josh coords
                            (SExcat1.ywin_image-0.5)-p[1], $
                            radius )
     if arg_present(err) then begin
        var=apphot( suberrimage^2, $
                    (SExcat1.xwin_image-0.5)-p[0], $  ;;SEx to josh coords
                    (SExcat1.ywin_image-0.5)-p[1], $
                    radius )
        errflux[igal,*] = sqrt(var)
     endif
  endfor
  if ~keyword_set(silent) then print
  mags = -2.5*alog10(fluxes)+(*self.zeropoint)[iband]
  if arg_present(err) then err = 2.5/2.*alog10((1.+errflux/fluxes)/(1.-errflux/fluxes))
  return, mags
end


;;;;;;;;;;;;;;;;;;;;;;
;;  wgals
;;
Function Cluster3::wgals, radius=radius
  s=(*self.SExcat)[*,0]
  c=self->color(ap_radius=10)
  wgals=where( s.flux_radius gt 2.2 $
               and s.mag_auto gt 17 and s.mag_auto lt 27 $
               and c gt -1 and c lt 2 $
               and s.a_image*s.kron_radius lt 200 )
  if n_elements(radius) ne 0 then begin
     ra = (*self.SExcat)[wgals,0].alphawin_j2000
     dec = (*self.SExcat)[wgals,0].deltawin_j2000
     ra0 = self.ra
     dec0 = self.dec
     r2 = (ra0-ra)^2*(cos(dec0*!dpi/180.d))^2 + (dec0-dec)^2
     r = sqrt(r2)*60. ;convert degree->arcmin
     wgals=wgals[where(r le radius)]
  endif
  return, wgals
end


;;;;;;;;;;;;;;;;;;;;;;
;;  wstars
;;
Function Cluster3::wstars, crange=crange
  s=(*self.SExcat)[*,0]
  wstars=where( s.flux_radius gt 1.4 and s.flux_radius lt 1.72 $
                and s.mag_auto gt 18 and s.mag_auto lt 24 )
  if n_elements(crange) ne 0 then begin
     c = self->color( ap_radius=10 )
     wstars=where( s.flux_radius gt 1.4 and s.flux_radius lt 1.72 $
                   and s.mag_auto gt 18 and s.mag_auto lt 24 $
                   and c ge crange[0] and c le crange[1] )
  endif
  return, wstars
end


;;;;;;;;;;;;;;;;;;;;;;
;;  IsophotSky
Function Cluster3::IsophotSky, gals=gals, $
                               band=band, $
                               store=store
  if ~ptr_valid(self.SExcat) || ~ptr_valid(self.galapagoscat) then begin
     message, $
        'IsophotSky requires a stored SExtractor catalog and stored GALAPAGOS catalog', $
        /continue
     return, 0
  endif
  if n_elements(gals) eq 0 then gals=indgen(n_elements(galapagoscat))
  store=keyword_set(store)

  if ~ptr_valid(self.skycat) then begin
     skycat = { skycat, $
                clusterid:   self.clusterid, $
                galid:       0L, $
                band:        '', $
                sky_value:   !values.d_nan, $
                sky_sigma:   !values.d_nan, $
                skyim_value: !values.d_nan, $
                exptime:     !values.d_nan, $
                nexp:        0 }
     sz=size(*self.SExcat, /dim)
     nobjs = sz[0]
     nbands = sz[1]
     skycat = replicate( skycat, nobjs, nbands )
     skycat.galid = rebin( lindgen(nobjs), nobjs, nbands )
     for iband=0, nbands-1 do skycat[*,iband].band = (*self.band)[iband]
     self.skycat = ptr_new(skycat)
  endif

  skycat = (*self.skycat)
  if n_elements(band) eq 0 then begin
     for iband=0, n_elements(*self.band)-1 do $
        skycat[gals,iband] = (self->IsophotSky(gals=gals, $
                                               band=(*self.band)[iband], $
                                               store=store ))[*,iband]
     return, skycat[gals,*]
  endif

  xband = (where(band eq (*self.band)))[0]
  for igal=0, n_elements(gals)-1 do begin
     counter, igal+1, n_elements(gals), 'IsophotSky on band '+band+': '
     gal1=(*self.galapagoscat)[gals[igal]]
     sex1=(*self.sexcat)[gals[igal],0]
     p = [gal1.xmin, gal1.ymin, gal1.xmax, gal1.ymax] + 250*[-1, -1, 1, 1]
     sciim = self->ImageSection( p, band=band, type='SCI' )
     mskim = self->ImageSection( p, band=band, type='MSK' )
     skyim = self->ImageSection( p, band=band, type='SKY' )
     whtim = self->ImageSection( p, band=band, type='WHT' )
     nxpim = self->ImageSection( p, band=band, type='NXP' )
     ellipse = { ellipse, $
                 center:[sex1.xwin_image-1, sex1.ywin_image-1]-p[0:1], $
                 majaxis:sex1.a_image*sex1.kron_radius, $
                 minaxis:sex1.b_image*sex1.kron_radius, $
                 theta:sex1.theta_image }
     w=where(~finite(mskim))
     skystruct = newisophotsky( sciim, mskim, ellipse, platescale=self.platescale )
     skycat[gals[igal],xband].sky_value=skystruct.sky_value
     skycat[gals[igal],xband].sky_sigma=skystruct.sky_sigma

     dist_ellipse, mask, [p[2]-p[0]+1, p[3]-p[1]+1], $
                   ellipse.center[0], ellipse.center[1], $
                   ellipse.majaxis/ellipse.minaxis, ellipse.theta-90.
     e=mask le ellipse.majaxis and finite(skyim)
     if total(e) eq 0 then begin
        skycat[gals[igal], xband].skyim_value = $
           skyim[ellipse.center[0], ellipse.center[1]]
        skycat[gals[igal], xband].exptime = $
           whtim[ellipse.center[0], ellipse.center[1]]
        skycat[gals[igal], xband].nexp = $
           nxpim[ellipse.center[0], ellipse.center[1]]
     endif else begin
        skycat[gals[igal], xband].skyim_value = biweight_mean(skyim[where(e ne 0.)])
        skycat[gals[igal], xband].exptime = biweight_mean(whtim[where(e ne 0.)])
        skycat[gals[igal], xband].nexp = round(biweight_mean(nxpim[where(e ne 0.)]))
     endelse
  endfor
  print
  if store then begin
     ptr_free, self.skycat
     self.skycat=ptr_new(skycat)
  endif
  return, skycat[gals, *]
end

;;;;;;;;;;;;;;;;;;;;;;
;;  SkyReg
Function Cluster3::SkyReg, gal, band
  if n_elements(band) eq 0 then band = 'z'

  gal1=(*self.galapagoscat)[gal]
  sex1=(*self.sexcat)[gal,0]
  p = [gal1.xmin, gal1.ymin, gal1.xmax, gal1.ymax] + 250*[-1, -1, 1, 1]
  sciim = self->ImageSection( p, band=band, type='SCI' )
  mskim = self->ImageSection( p, band=band, type='MSK' )
  ellipse = { ellipse, $
              center:[sex1.xwin_image-1, sex1.ywin_image-1]-p[0:1], $
              majaxis:sex1.a_image*sex1.kron_radius, $
              minaxis:sex1.b_image*sex1.kron_radius, $
              theta:sex1.theta_image }
  junk = newisophotsky( sciim, mskim, ellipse, platescale=self.platescale, skyreg=skyreg )
  return, skyreg
end


;;;;;;;;;;;;;;;;;;;;;;
;;  GetStars
Function Cluster3::GetStars, band=band, $
                             starsize=starsize, $
                             rthresh=rthresh, $
                             minstars=minstars, $
                             _extra=extra
  if n_elements(band) eq 0 then band=(*self.band)[0]
  xband = (where(band eq *self.band))[0]
  if n_elements(starsize) eq 0 then begin
     starsize=fix(41)
     if fix(starsize/2.) eq starsize/2. then starsize += 1 ;;ensure odd
  endif
  innersize = fix(starsize*0.23)
  if fix(innersize/2.) eq innersize/2. then innersize += 1 ;;ensure odd

  wstars = self->wstars( _extra=extra )

  if n_elements(rthresh) ne 0 then begin
     ;; radial threshold in arcmin.
     ra = (*self.SExcat)[wstars,0].alphawin_j2000
     dec = (*self.SExcat)[wstars,0].deltawin_j2000
     ra0 = self.ra
     dec0 = self.dec
     r2 = (ra0-ra)^2*(cos(dec0*!dpi/180.d))^2 + (dec0-dec)^2
     r = sqrt(r2)*60. ;;convert degree->arcmin
     if n_elements(minstars) eq 0 then begin
        wstars=wstars[where(r le rthresh)]
     endif else begin
        if n_elements(wstars) lt minstars then begin
           print, 'not enough stars!'
           print, 'using all stars found!'
        endif else begin
           if total(r le rthresh) lt minstars then begin
              s=sort(r)
              wstars=wstars[s[0:minstars-1]]
           endif else begin
              wstars=wstars[where(r le rthresh)]
           endelse
        endelse
     endelse
  endif

  halfsize=(starsize-1)/2
  onestar = { star, $
              clusterid:   '', $
              galid:       0L, $
              image:       fltarr(starsize,starsize), $
              model:       fltarr(innersize,innersize), $
              errim:       fltarr(starsize,starsize), $
              modelparams: fltarr(7), $
              x_SEx:       0.d, $
              y_SEx:       0.d }
  for iw=0, n_elements(wstars)-1 do begin
     star = onestar
     star.clusterid=self.clusterid
     star.galid=wstars[iw]
     p = round([(*self.SExcat)[wstars[iw],0].xwin_image-1, $ ;; SEx coords -> IDL coords
                (*self.SExcat)[wstars[iw],0].ywin_image-1, $
                (*self.SExcat)[wstars[iw],0].xwin_image-1, $
                (*self.SExcat)[wstars[iw],0].ywin_image-1])+[-1,-1,1,1]*halfsize
     star.image = self->ImageSection( p, band=band, type='SCI' )
     star.errim = self->ImageSection( p, band=band, type='ERR' )
     if n_elements(star.errim) eq 0 || total(~finite(star.errim)) ne 0 then continue
     sky_value = (*self.skycat)[wstars[iw],xband].sky_value
     if finite(sky_value) then star.image -= sky_value
     star.x_SEx = (*self.SExcat)[wstars[iw],0].xwin_image $
                  - fix((*self.SExcat)[wstars[iw],0].xwin_image)
     star.y_SEx = (*self.SExcat)[wstars[iw],0].ywin_image $
                  - fix((*self.SExcat)[wstars[iw],0].ywin_image)
     lowerlimit = (starsize-innersize)/2
     upperlimit = lowerlimit+innersize-1
     inner = (star.image)[lowerlimit:upperlimit,lowerlimit:upperlimit]
     innererr = (star.errim)[lowerlimit:upperlimit,lowerlimit:upperlimit]
     xs=findgen(innersize)+lowerlimit+0.5
     ys=findgen(innersize)+lowerlimit+0.5
     star.model = mpfit2dpeak(inner, a, xs, ys, error=innererr)
     star.modelparams = a
     if n_elements(stars) eq 0 then stars=star else stars=[stars,star]
  endfor
  if n_elements(stars) eq 0 then return, -1
  return, stars
end


;;;;;;;;;;;;;;;;;;;;;;
;;  Galfit
Function Cluster3::Galfit, gals=gals, $
                           outdir=outdir, $
                           take2=take2, $
                           noerr=noerr, $
                           doall=doall, $
                           store=store, $
                           _extra=extra
  if ~ptr_valid(self.galapagoscat) || ~ptr_valid(self.SExcat) $
     || ~ptr_valid(self.skycat) || ~ptr_valid(self.sciim) $
     || ~ptr_valid((*self.sciim)[0]) || ~ptr_valid(self.psf) $
     || ~ptr_valid((*self.psf)[0]) then begin
     message, 'Galfit requires a stored SExtractor catalog, a stored GALAPAGOS catalog,' $
              + ' a stored sky background catalog, and stored images/PSFs', /continue
     return, 0
  endif

  if n_elements(outdir) eq 0 then outdir=self.rootdir+'galfit/'
  outdir=directoryify(outdir)
  if n_elements(gals) eq 0 then gals = lindgen( n_elements( *self.galapagoscat ) )
  if total(gals) eq -1 then return, 0
  take2 = keyword_set(take2)
  noerr = keyword_set(noerr)
  store = keyword_set(store)
  doall = keyword_set(doall)

  galfitcat1 = { galfitcat, $
                 clusterid:  self.clusterid, $
                 galid:      0L, $
                 xpos:       !values.d_nan, $
                 ypos:       !values.d_nan, $
                 mag:        !values.f_nan, $
                 Re:         !values.f_nan, $
                 boa:        !values.f_nan, $
                 pa:         !values.f_nan, $
                 n:          !values.f_nan, $
                 errxpos:    !values.f_nan, $
                 errypos:    !values.f_nan, $
                 errmag:     !values.f_nan, $
                 errRe:      !values.f_nan, $
                 errboa:     !values.f_nan, $
                 errpa:      !values.f_nan, $
                 errn:       !values.f_nan, $
                 chi2:       !values.f_nan, $
                 dof:        -1, $
                 chi2perdof: !values.f_nan, $
                 galfittime: !values.f_nan, $
                 galfitretry: -1 }
  
  if ~ptr_valid(self.galfitcat) then begin
     galfitcat = replicate( galfitcat1, n_elements(*self.galapagoscat) )
     galfitcat.galid = (*self.galapagoscat).galid
     self.galfitcat = ptr_new(galfitcat)
  endif


  ;; still not sure if this works!!!
  if doall then begin
     (*self.galfitcat)[gals] = replicate( galfitcat1, n_elements(gals))
     (*self.galfitcat)[gals].galid = gals
     galfitcat0 = self->Galfit( gals=gals, outdir=outdir, store=store, _extra=extra )
     wfail = where( finite(galfitcat0.galfittime) and ~finite(galfitcat0.re) )
     if wfail[0] ne -1 then begin
        galfitcat1 = self->Galfit( gals=gals[wfail], outdir=outdir, $
                                   /take2, store=store, _extra=extra )
        wfail2 = where( finite(galfitcat1.galfittime) and ~finite(galfitcat1.re) )
        if wfail2[0] ne -1 then begin
           galfitcat2 = self->Galfit( gals=gals[wfail[wfail2]], outdir=outdir, $
                                      /noerr, store=store, _extra=extra )
           wfail3 = where( finite(galfitcat2.galfittime) and ~finite(galfitcat2.re) )
           if wfail3[0] ne -1 then begin
              galfitcat3 = self->Galfit( gals=gals[wfail[wfail2[wfail3]]], outdir=outdir, $
                                         /take2, /noerr, store=store, _extra=extra )
           endif
        endif
     endif
     if wfail[0] ne -1 then galfitcat0[wfail]=galfitcat1
     if n_elements(wfail2) ne 0 && wfail2[0] ne -1 then galfitcat0[wfail[wfail2]]=galfitcat2
     if n_elements(wfail3) ne 0 && wfail3[0] ne -1 then galfitcat0[wfail[wfail2[wfail3]]]=galfitcat3
     return, galfitcat0
  endif else begin

     galfitcat = (*self.galfitcat)
     case 1 of
        (~keyword_set(take2) and ~keyword_set(noerr)) : $
           print, 'Galfit attempt 1'
        (keyword_set(take2) and ~keyword_set(noerr))  : $
           print, 'Galfit attempt 2'
        (~keyword_set(take2) and keyword_set(noerr))  : $
           print, 'Galfit attempt 3'
        (keyword_set(take2) and keyword_set(noerr))   : $
           print, 'Galfit attempt 4'
     endcase 
     for igal=0L, n_elements(gals)-1 do begin
        counter, igal+1, n_elements(gals), 'Running galfit: '
        glp1 = (*self.galapagoscat)[gals[igal]]
        starttime = systime(1)

        galfit1 = self->GalfitObj( gals[igal], band=band, outdir=outdir, $
                                   take2=take2, noerr=noerr, _extra=extra )

        if size(galfit1, /tname) ne 'OBJREF' then begin
           galfitcat[gals[igal]].galfittime = systime(1)-starttime
           continue
        endif
        galfit1->Execute
        results = galfit1->Results()
        if size(results, /tname) eq 'STRUCT' then begin
           results.xpos += glp1.xmin
           results.ypos += glp1.ymin
           dest = galfitcat[gals[igal]]
           struct_assign, results, dest, /nozero
           galfitcat[gals[igal]] = dest
           case 1 of 
              (~keyword_set(take2) and ~keyword_set(noerr)) : $
                 galfitcat[gals[igal]].galfitretry=0
              (keyword_set(take2) and ~keyword_set(noerr))  : $
                 galfitcat[gals[igal]].galfitretry=1
              (~keyword_set(take2) and keyword_set(noerr))  : $
                 galfitcat[gals[igal]].galfitretry=2
              (keyword_set(take2) and keyword_set(noerr))   : $
                 galfitcat[gals[igal]].galfitretry=3
           endcase        
        endif
        galfitcat[gals[igal]].galfittime = systime(1)-starttime
        obj_destroy, galfit1
     endfor
     print
     if keyword_set(store) then begin
        ptr_free, self.galfitcat
        self.galfitcat = ptr_new(galfitcat)
     endif
  endelse
  return, galfitcat[gals]
end


;;;;;;;;;;;;;;;;;;;;;;
;;  GalfitObj
Function Cluster3::GalfitObj, gal, $
                              band=band, $
                              outdir=outdir, $
                              _extra=extra
;; extra includes:
;; noerr=noerr
;; take2=take2 for 2nd chance GALFITing.
;; sn = {x, y, mag} for PSF fitting

  if n_elements(band) eq 0 then band=(*self.band)[0]
  xband = (where(band eq (*self.band)))[0]
  if n_elements(outdir) eq 0 then outdir=rootdir+'galfit/'
  outdir=directoryify(outdir)

  psf = float(undersampleimage(*(*self.psf)[xband], self.samplefactor))
  psf /= max(psf)

  hdr = { exptime:   (*self.skycat)[gal,xband].exptime, $
          gain:      1.0, $
          rdnoise:   5.0, $
          ncombine:  (*self.skycat)[gal,xband].nexp, $
          zeropoint: (*self.zeropoint)[xband], $
          band:band }

  sciim = self->galimage( gal, band=band, type='SCI' )
  errim = self->galimage( gal, band=band, type='ERR' )

  galfit1 = MakeGalfit3( sciim, errim, psf, $
                         (*self.sexcat)[*,xband], $
                         *self.galapagoscat, (*self.skycat)[*,xband], $
                         gal, hdr, outdir=outdir, _extra=extra )

  return, galfit1
end


;;;;;;;;;;;;;;;;;;;;;;
;;  RecallGalfit
Function Cluster3::RecallGalfit, $
   gal, $
   outdir=outdir

  if n_elements(outdir) eq 0 then outdir=self.rootdir+'galfit/'
  outdir=directoryify(outdir)
  prefix=self.clusterid+string(gal,format='(I05)')+'z'
  filename=outdir+prefix+'.fits'
  if ~(file_info(filename)).exists then return, -1.
  struct = mrdfits( filename, 1, /silent )
  return, struct
end


;;;;;;;;;;;;;;;;;;;;;;
;;  Morphology
Function Cluster3::Morphology, $
   gals=gals, $
   store=store, $
   segtype=segtype, $
   skyfactor=skyfactor, $
   rotresid=rotresid, $
   new_mask=new_mask, $
   imagetype=imagetype, $
   _extra=extra

  if n_elements(segtype) eq 0 then segtype=1 ;;Josh segmentation
  if n_elements(imagetype) eq 0 then imagetype=2
  if n_elements(skyfactor) eq 0 then skyfactor=1.5
  store=keyword_set(store)
  morphcat1 = { morphcat, $
                clusterid:  self.clusterid, $
                galid:      0L, $
                qpflux:     !values.f_nan, $
                asym:       !values.f_nan, $
                asymerr:    !values.f_nan, $
                gini:       !values.f_nan, $
                ginierr:    !values.f_nan }
  
  if ~ptr_valid(self.morphcat) then begin
     morphcat = replicate( morphcat1, n_elements(*self.galapagoscat) )
     morphcat.galid = (*self.galapagoscat).galid
     self.morphcat = ptr_new(morphcat)
  endif
  
  morphcat=*self.morphcat
  
  for igal=0, n_elements(gals)-1 do begin
     counter, igal+1, n_elements(gals), 'Measuring morphology: '
     case imagetype of
        1: begin
           image = self->galimage(gals[igal])-(*self.skycat)[gals[igal],0].sky_value
        end
        2: begin
           struct = self->recallgalfit( gals[igal], _extra=extra )
           if size(struct, /tname) ne 'STRUCT' then image=self->galimage(gals[igal])*0. $
           else image = struct.intsb
        end
     endcase
     if total(image ne 0) eq 0 then continue
     skysig = (*self.skycat)[gals[igal],0].sky_sigma
     sz=size( image, /dim )
     center=sz/2
     case segtype of 
        1: begin ;;josh style segmentation
           seg = ACSsegmask( image, skyfactor*skysig )
           mask0 = seg eq seg[center[0], center[1]]
        end
        2: begin ;;hot SExtractor segmentation
           seg = self->galimage( gals[igal], type='SEGHOT' )
           mask0 = seg eq seg[center[0], center[1]]
        end
        3: begin ;;hot/cold SExtractor segmentation
           if (*self.sexcat)[gals[igal],0].cold then seg = self->galimage( gals[igal], type='SEGCLD' ) $
           else seg = self->galimage( gals[igal], type='SEGHOT' )
           mask0 = seg eq seg[center[0], center[1]]
        end
        4: begin ;;use the interloper subtracted image to segment
           struct = self->recallgalfit( gals[igal], _extra=extra )
           seg = ACSsegmask( struct.intsb, skyfactor*skysig )
           mask0 = seg eq seg[center[0], center[1]]
        end
     endcase 
     
     ;;grow mask 2 pixels
     mask = mask0
     for i=-2,2 do begin
        for j=-2,2 do begin
           if abs(i)+abs(j) gt 2 then continue
           mask = mask or shift( mask0, i, j )
        endfor
     endfor
     
     if arg_present(rotresid) then $
        morph = morph( image, mask, skysig, /symmetric, $
                       new_mask=new_mask, rotresid=rotresid, eta=0.25 ) $ ;;obtain morphology!
     else $
        morph = morph( image, mask, skysig, /symmetric, eta=0.25 )
     struct_assign, morph, morphcat1
     morphcat[gals[igal]] = morphcat1
     morphcat[gals[igal]].galid=gals[igal]
  endfor
  print

  if store then begin
     ptr_free, self.morphcat
     self.morphcat=ptr_new(morphcat)
  endif
  return, morphcat[gals]
end


;;;;;;;;;;;;;;;;;;;;;;
;;  Masks - just a diagnostic tool for Josh....
Function Cluster3::masks, gal
  skyfactor=1.5
  image = self->galimage(gal)-(*self.skycat)[gal,0].sky_value
  skysig = (*self.skycat)[gal,0].sky_sigma
  sz=size( image, /dim )
  center=sz/2
  seg = ACSsegmask( image, skyfactor*skysig )
  mask1 = seg eq seg[center[0], center[1]]
  seg = self->galimage( gal, type='SEGHOT' )
  mask2 = seg eq seg[center[0], center[1]]
  if (*self.sexcat)[gal,0].cold then seg = self->galimage( gal, type='SEGCLD' ) $
  else seg = self->galimage( gal, type='SEGHOT' )
  mask3 = seg eq seg[center[0], center[1]]
  return, [image, mask1, mask2, mask3]
end


;;;;;;;;;;;;;;;;;;;;;;
;;  XPSFCat
Function Cluster3::XPSFcat, $
   gals=gals, $
   band=band, $
   store=store

  if n_elements(gals) eq 0 then gals=lindgen(n_elements(*self.galfitcat))
  if total(gals) eq -1 then return, 0
  store=keyword_set(store)
  if ~ptr_valid(self.xPSFcat) then begin
     xPSFcat1 = { xPSFcat, $
                  clusterid:  self.clusterid, $
                  galid:      0L, $
                  band:       '', $
                  xpsfap:     fltarr(5)*!values.f_nan, $
                  xpsffac:    fltarr(5)*!values.f_nan }
     xPSFcat = replicate( xPSFcat1, n_elements(*self.galfitcat), n_elements(*self.band) )
     xPSFcat.galid = (*self.SExcat).galid
     xPSFcat.band = (*self.SExcat).band
     self.xPSFcat = ptr_new(xPSFcat)
  endif
  xPSFcat = *self.xPSFcat

  if n_elements(band) eq 0 then begin
     for iband=0, n_elements(*self.band)-1 do begin
        band=(*self.band)[iband]
        xPSFcat[gals,iband] = (self->xPSFcat(gals=gals, band=band, store=store ))[*,iband]
     endfor
     return, xPSFcat[gals,*]
  endif

  iband=(where(band eq (*self.band)))[0]
  for igal=0, n_elements(gals)-1 do begin
     counter, igal+1, n_elements(gals), 'Creating XPSFcat for band '+band+': '
     gal1 = (*self.galfitcat)[gals[igal]]
     SEx1 = (*self.SExcat)[gals[igal],0]
     sky1 = (*self.skycat)[gals[igal],iband]
     max_radius = 20. > 1.5*gal1.re < 500.
     
     p = floor([ SEx1.xwin_image-1, $ ;;SEx coords to IDL coords
                 SEx1.ywin_image-1, $
                 SEx1.xwin_image-1, $
                 SEx1.ywin_image-1 ]) + [-1,-1,1,1]*floor((max_radius*1.2) + 20)
     cvlim = self->ImageSection(p, band=band, type='CVL')
     if finite(sky1.sky_value) then cvlim -= sky1.sky_value
     apfluxes = apphot( cvlim, $
                        SEx1.xwin_image-0.5-p[0], $  ;;SEx to josh coords
                        SEx1.ywin_image-0.5-p[1], $
                        [3., 5., 10., 15., 20.] )
     facfluxes = apphot( cvlim, $
                         SEx1.xwin_image-0.5-p[0], $
                         SEx1.ywin_image-0.5-p[1], $
                         [0.5, 0.75, 1.0, 1.25, 1.5]*(gal1.re+1.2) ) ;;correction for PSF convol.
     xPSFcat[gals[igal],iband].xpsfap = -2.5*alog10(apfluxes)+(*self.zeropoint)[iband]
     xPSFcat[gals[igal],iband].xpsffac = -2.5*alog10(facfluxes)+(*self.zeropoint)[iband]
  endfor
  print
  if keyword_set(store) then begin
     ptr_free, self.xPSFcat
     self.xPSFcat = ptr_new(xPSFcat)
  endif
  return, xPSFcat[gals,*]
end


;;;;;;;;;;;;;;;;;;;;;;
;;  AddSpectroscopy
Pro Cluster3::AddSpectroscopy, $
   dir, $  ;; directory to add
   filename
  dir=directoryify(dir)
  spectra = addspectra1( dir, filename, (*self.SExcat)[*,0] )
  if ptr_valid(self.spectra) then begin
     oldspectra=*self.spectra
     spectra = [oldspectra, spectra]
     ptr_free, self.spectra
     self.spectra = ptr_new(spectra)
  endif else self.spectra = ptr_new(spectra)
  self->UpdateSpecCat
end


;;;;;;;;;;;;;;;;;;;;;;
;;  MeasureO2
Pro Cluster3::MeasureO2, plot=plot1
  resolve_routine, 'spectrum__define', /no_recompile
  resolve_routine, 'dopplershift__define', /no_recompile
  id=self.clusterid
  if ~ptr_valid(self.spectra) then return
  spectra=*self.spectra
  s=sort(spectra.galid)
  u=uniq(spectra.galid, s)
  for i=0, n_elements(u)-1 do begin
     w=where(spectra.galid eq spectra[u[i]].galid)
     for iw=0, n_elements(w)-1 do begin
        counter, iw+1, n_elements(w), 'Object '+str(i+1)+' of '+str(n_elements(u))+'. Spec '
        s0=spectra[w[iw]]
        s1=s0.spectrum
        galid=s0.galid
        num=string(galid, format='(I05)')
        if iw lt 26 then name=id+num+string(byte(97+iw)) $
        else name=id+num+string(byte(97+(iw/26)))+string(byte(97+(iw mod 26)))
        if obj_valid(s1) then begin
           if keyword_set(plot1) then begin
              plot='/home/scpdata02/cluster3/diagnostic/o2/'+name+'.o2.eps'
              measureo2, s1, O2EWval=O2_EW, O2EWerr=O2_EWerr, plot=plot
           endif else measureo2, s1, O2EWval=O2_EW, O2EWerr=O2_EWerr
           spectra[w[iw]].o2_ew=O2_EW
           spectra[w[iw]].o2_ewerr=O2_EWerr
        endif
     endfor
  endfor
  ptr_free, self.spectra
  self.spectra=ptr_new(spectra)
  print
end


;;;;;;;;;;;;;;;;;;;;;;
;;  YoII
Pro Cluster3::YoII
  readcol, '/home/jmeyers314/scp1/CL-spec/Literature/DeMarco2007/Ydata.dat', $
           format='A,A,A,D,D,X,X,D,D,D,D,D,D', $
           id, ra0, dec0, z, zerr, o2flux, o2fluxerr, o2ew, o2ewerr, sfr, sfrerr
  spectra=*self.spectra
  for i=0, n_elements(ra0)-1 do begin
     radeg=ten(double(strsplit( ra0[i], ':', /extract )))*15.+3.815e-5
     decdeg=ten(double(strsplit( dec0[i], ':', /extract )))+3.898e-6
     dist2=(radeg-spectra.ra)^2*cos(decdeg*!pi/180.)^2+(decdeg-spectra.dec)^2
     dist2 *= 3600.d^2
     dist = sqrt(dist2)
     m=min( dist, w )

     print, m
     print, radeg, spectra[w].ra, radeg-spectra[w].ra
     print, decdeg, spectra[w].dec, decdeg-spectra[w].dec
     print
     print
     spectra[w].o2_ew = o2ew[i]
     spectra[w].o2_ewerr = o2ewerr[i]
  endfor
  ptr_free, self.spectra
  self.spectra=ptr_new(spectra)
  
end


;;;;;;;;;;;;;;;;;;;;;;
;;  UpdateSpecCat
;;  Takes spectra catalog and uses it to fill in galaxy [OII], z,
;;  zerr, zqual, comment, etc.
Pro Cluster3::UpdateSpecCat
  speccat = { SpecCat, $
              galid:    0L, $
              z:        -2.d, $
              zerr:     -1.d, $
              zqual:    'N', $
              comment:  '', $
              o2_EW:    !values.f_nan, $
              o2_EWerr: !values.f_nan }
  speccat = replicate( speccat, n_elements(*self.galfitcat) )
  speccat.galid=(*self.galfitcat).galid
  if ptr_valid(self.spectra) then begin ;;find galid's for existing spectra
     gals = ((*self.spectra).galid)[uniq( (*self.spectra).galid, $
                                          sort( (*self.spectra).galid ) )]
  endif
  for igal=0, n_elements(gals)-1 do begin
     if gals[igal] lt 0 then continue;; these are the non-ACS spectra...
     w = where( (*self.spectra).galid eq gals[igal] )
     wgood = where( (*self.spectra)[w].z ne -1. )
     if wgood[0] eq -1 then begin  ;;didn't get any redshifts for this object
        speccat[gals[igal]].z=-1.d
        speccat[gals[igal]].zerr=-1.d
        speccat[gals[igal]].zqual='F'
     endif else begin  ;;got at least 1 redshift...
        ;; take average z
        speccat[gals[igal]].z=mean( (*self.spectra)[w[wgood]].z ) 
        ;;take minimum zerr
        speccat[gals[igal]].zerr=min( (*self.spectra)[w[wgood]].zerr ) 
        ;;take best zqual
        speccat[gals[igal]].zqual=(((*self.spectra)[w[wgood]].zqual) $
                                   [sort( (*self.spectra)[w[wgood]].zqual )])[0] 
        ;;take o2_EW from best zqual
        speccat[gals[igal]].o2_EW=(((*self.spectra)[w[wgood]].o2_EW) $
                                   [sort( (*self.spectra)[w[wgood]].zqual )])[0] 
        speccat[gals[igal]].o2_EWerr=(((*self.spectra)[w[wgood]].o2_EWerr) $
                                      [sort( (*self.spectra)[w[wgood]].zqual )])[0] 
        ;; unionize comments
        comment = strjoin((*self.spectra)[w[wgood]].comment+',') 
        comment = strsplit(comment, ',', /extract)
        speccat[gals[igal]].comment=strjoin( comment[uniq( comment, $
                                                           sort(comment) )],',' )
     endelse
  endfor
  ptr_free, self.speccat
  self.speccat = ptr_new(speccat)
  ;; organize the non-ACS spectra
  if ~ptr_valid(self.spectra) then return
  spectra = *self.spectra
  w=where(spectra.galid lt 0)
  if w[0] eq -1 then return
  completed = spectra.galid ge 0
  curid = -1
  while total(completed) ne n_elements(spectra) do begin
     notcompleted = where(completed eq 0)
     nextgal = notcompleted[0]
     ra0 = spectra[nextgal].ra
     dec0 = spectra[nextgal].dec
     ra = spectra[notcompleted].ra
     dec = spectra[notcompleted].dec
     dist2=(ra-ra0)^2*(cos(dec*!dpi/180.d))^2+(dec-dec0)^2
     near = where(dist2 le 1./3600./3600.)
     spectra[notcompleted[near]].galid = curid
     curid -= 1
     completed[notcompleted[near]]=1
  endwhile
  ptr_free, self.spectra
  self.spectra = ptr_new(spectra)
end


Function Cluster3::CompleteSpecCat
  speccat = { CompleteSpecCat, $
              galid:    0L, $
              ra:       !values.d_nan, $
              dec:      !values.d_nan, $
              z:        -2.d, $
              zerr:     -1.d, $
              zqual:    'N', $
              comment:  '', $
              o2_EW:    !values.f_nan, $
              o2_EWerr: !values.f_nan }
  if ~ptr_valid(self.spectra) then return, 0.
  spectra = *self.spectra
  gals = (spectra.galid)[uniq( spectra.galid, $
                               sort( spectra.galid ) )]

  speccat = replicate( speccat, n_elements(gals) )
  speccat.galid = gals
  for igal=0, n_elements(gals)-1 do begin
     w = where( spectra.galid eq gals[igal] )
     wgood = where( spectra[w].z ne -1. )
     if wgood[0] eq -1 then begin  ;;didn't get any redshifts of this object
        speccat[igal].z=-1.d
        speccat[igal].zerr=-1.d
        speccat[igal].zqual='F'
     endif else begin  ;;got at least 1 redshift...
        speccat[igal].ra=mean( spectra[w[wgood]].ra )
        speccat[igal].dec=mean( spectra[w[wgood]].dec )
        ;; take average z
        speccat[igal].z=mean( spectra[w[wgood]].z ) 
        ;;take minimum zerr
        speccat[igal].zerr=min( spectra[w[wgood]].zerr ) 
        ;;take best zqual
        speccat[igal].zqual=((spectra[w[wgood]].zqual) $
                             [sort( spectra[w[wgood]].zqual )])[0] 
        ;;take o2_EW from best zqual
        speccat[igal].o2_EW=((spectra[w[wgood]].o2_EW) $
                             [sort( spectra[w[wgood]].zqual )])[0] 
        speccat[igal].o2_EWerr=((spectra[w[wgood]].o2_EWerr) $
                                [sort( spectra[w[wgood]].zqual )])[0] 
        ;; unionize comments
        comment = strjoin(spectra[w[wgood]].comment+',') 
        comment = strsplit(comment, ',', /extract)
        speccat[igal].comment=strjoin( comment[uniq( comment, $
                                                     sort(comment) )],',' )
     endelse
  endfor
  return, speccat
end

;;;;;;;;;;;;;;;;;;;;;;
;;  Summary
Function Cluster3::Summary
  commontags = ['CLUSTERID','GALID','XWIN_IMAGE','YWIN_IMAGE', $
                'ALPHAWIN_J2000','DELTAWIN_J2000', $
                'A_IMAGE','B_IMAGE','THETA_IMAGE', 'ELLIPTICITY', $
                'KRON_RADIUS','PETRO_RADIUS','COLD']
  outcat = struct_trimtags((*self.SExcat)[*,0], select=commontags)
  nbands = n_elements(*self.band)
  for iband=0, nbands-1 do begin
     SExcat = struct_trimtags((*self.SExcat)[*,iband], except=[commontags,'BAND'])
     tags = tag_names(SExcat)
     SExcat = struct_renametags(SExcat, (*self.band)[iband]+tags)
     outcat = struct_addtags(outcat, SExcat)
  endfor
  if ptr_valid(self.galapagoscat) then begin
     galapagoscat = struct_trimtags(*self.galapagoscat, except=['CLUSTERID','GALID'])
     outcat = struct_addtags(outcat, galapagoscat)
  endif
  if ptr_valid(self.skycat) then begin
     for iband=0, nbands-1 do begin
        skycat = struct_trimtags((*self.skycat)[*,iband], $
                                 except=['CLUSTERID','GALID','BAND'])
        tags = tag_names(skycat)
        skycat = struct_renametags(skycat, (*self.band)[iband]+tags)
        outcat = struct_addtags(outcat, skycat)
     endfor
  endif
  if ptr_valid(self.galfitcat) then begin
     galfitcat = struct_trimtags(*self.galfitcat, except=['CLUSTERID','GALID'])
     outcat = struct_addtags(outcat, galfitcat)
  endif
  if ptr_valid(self.morphcat) then begin
     morphcat = struct_trimtags(*self.morphcat, except=['CLUSTERID','GALID'])
     outcat = struct_addtags(outcat, morphcat)
  endif
;  if ptr_valid(self.morph2cat) then begin
;     morph2cat = struct_trimtags(*self.morph2cat, except=['CLUSTERID','GALID'])
;     outcat = struct_addtags(outcat, morph2cat)
;  endif
  if ptr_valid(self.CLEANcat) then begin
     for iband=0, nbands-1 do begin
        CLEANcat = struct_trimtags((*self.CLEANcat)[*,iband], $
                                   except=['CLUSTERID','GALID','BAND'])
        tags = tag_names(CLEANcat)
        CLEANcat = struct_renametags(CLEANcat, (*self.band)[iband]+tags)
        outcat = struct_addtags(outcat, CLEANcat)
     endfor
  endif
  if ptr_valid(self.xPSFcat) then begin
     for iband=0, nbands-1 do begin
        xPSFcat = struct_trimtags((*self.xPSFcat)[*,iband], $
                                  except=['CLUSTERID','GALID','BAND'])
        tags = tag_names(xPSFcat)
        xPSFcat = struct_renametags(xPSFcat, (*self.band)[iband]+tags)
        outcat = struct_addtags(outcat, xPSFcat)
     endfor
  endif
  if ptr_valid(self.speccat) then begin
     speccat = struct_trimtags(*self.speccat, except=['CLUSTERID','GALID'])
     outcat = struct_addtags(outcat, speccat)
  endif
  return, outcat
end


;;;;;;;;;;;;;;;;;;;;;;
;;  CreateRegion

Pro Cluster3::CreateRegion, filename, gals=gals
  if n_elements(gals) eq 0 then gals=lindgen( ( size( *self.SExcat, /dim ) )[0] )
  openw, lun, filename, /get_lun
  printf, lun, 'global color=green font="helvetica 10 normal" '+ $
          'select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source'
  printf, lun, 'fk5'
  printf, lun
  printf, lun
  for i=0, n_elements(gals)-1 do begin
     if (*self.SExcat)[gals[i],0].cold eq 1 then color='blue' $
     else color = 'red'
     scat = (*self.SExcat)[gals[i],0]
     printf, lun, 'ellipse(' $
             +string(scat.alphawin_j2000,format='(D11.5)')+',' $
             +string(scat.deltawin_j2000,format='(D11.5)')+',' $
             +string(scat.a_image*scat.kron_radius*0.05,format='(D8.2)')+'",' $
             +string(scat.b_image*scat.kron_radius*0.05,format='(D8.2)') $
             +'", '+string(scat.theta_image,format='(D8.2)') $
             +') # tag={aperture} text={' $
             +string(scat.galid,format='(I05)')+'} color='+color
     printf, lun
  endfor
  close, lun
  free_lun, lun
end

;;;;;;;;;;;;;;;;;;;;;;
;;  NearestObjs
Function Cluster3::NearestObjs, ra, dec, dist=dist
  SExcat=(*self.SExcat)[*,0]
  dist2=(ra-SExcat.alphawin_j2000)^2*(cos(dec*!dpi/180.d))^2+(dec-SExcat.deltawin_j2000)^2
  s=sort(dist2)
  if arg_present(dist) then begin
     dist = sqrt(dist2[s[0:9]])*3600.d ;; distances in arcsec.
  endif
  return, SExcat[s[0:9]].galid
end


;;;;;;;;;;;;;;;;;;;;;;
;;  SNhost
Function Cluster3::SNhost, snname, dist=dist
  if (file_info(self.rootdir+'../supernova.txt')).exists then begin
     readcol, self.rootdir+'../supernova.txt', $
              snnames, snnick, z, type, ra, dec, hostra, hostdec, $
              format='A,A,A,A,A,A,A,A', /silent
     w=where(snname eq snnames or snname eq snnick)
     if w[0] eq -1 then return, -1
     get_coords, coords, instring=hostra[w[0]]+' '+hostdec[w[0]]
  endif else begin
     readcol, self.rootdir+'../GOODSsne.txt', $
              name, nick, ra, dec, format='A,A,A,A,X', /silent
     w=where(snname eq nick or snname eq name)
     if w[0] eq -1 then return, -1
     get_coords, coords, instring=ra[w[0]]+' '+dec[w[0]]
  endelse
  ra_deg=coords[0]*15.
  dec_deg=coords[1]
  if arg_present(dist) then begin ;; distance to host (nearest catalog object)
     out = (self->nearestobjs( ra_deg, dec_deg, dist=dist ))[0]
     dist = dist[0]
     return, out
  endif
  return, (self->nearestobjs( ra_deg, dec_deg ))[0]
end


;;;;;;;;;;;;;;;;;;;;;;
;;  SNcoord
Function Cluster3::SNcoord, snname, name1=name1
  if (file_info(self.rootdir+'../supernova.txt')).exists then begin
     readcol, self.rootdir+'../supernova.txt', $
              snnames, snnick, z, type, ra, dec, hostra, hostdec, $
              format='A,A,A,A,A,A,A,A', /silent
     w=where(snname eq snnames or snname eq snnick)
     if w[0] ne -1 then name1=snnames[w[0]]
  endif else begin
     readcol, self.rootdir+'../GOODSsne.txt', $
              name, nick, ra, dec, z, format='A,A,A,A,A'
     w=where( snname eq name or snname eq nick )
     if w[0] ne -1 then name1=name[w[0]]
  endelse
  if w[0] eq -1 then return, -1
  get_coords, coords, instring=ra[w[0]]+' '+dec[w[0]]
  ra_deg=coords[0]*15.
  dec_deg=coords[1]
  hdr = *self.hdr
  adxy, hdr, ra_deg, dec_deg, xpix, ypix
  return, [xpix, ypix]
end


;;;;;;;;;;;;;;;;;;;;;;
;;  VisMorphStruct
Function Cluster3::VisMorphStruct,  $
   filename, $
   radec=radec

  VisMorphCat0 = { VisMorphCat, $
                   clusterid: self.clusterid, $
                   galid:     0L, $
                   morph:     '', $
                   ttype:     99, $
                   gini:      !values.f_nan, $
                   ginierr:   !values.f_nan, $
                   asym:      !values.f_nan, $
                   asymerr:   !values.f_nan, $
                   color:     !values.f_nan }
  VisMorphCat = replicate( VisMorphCat0, n_elements(*self.galfitcat) )
  VisMorphCat.galid = (*self.galfitcat).galid
  VisMorphCat.gini = (*self.morphcat).gini
  VisMorphCat.ginierr = (*self.morphcat).ginierr
  VisMorphCat.asym = (*self.morphcat).asym
  VisMorphCat.asymerr = (*self.morphcat).asymerr
  VisMorphCat.color = self->color(/xpsf)

  if ~keyword_set(radec) then begin
     if self.clusterid eq 'F' then offset = [-1144.5, -1069]
     if self.clusterid eq 'H' then offset = [-1013.5, -923]
     fmt='I4,I6,F9.2,F9.2,F6.2,F7.2,I3,A100'
     readfmt, filename, fmt, id, id2, x, y, mag, blah, t, morf
     x += offset[0]
     y += offset[1]
     s=(*self.SExcat)[*,0]
     for iid=0, n_elements(id)-1 do begin
        dist2 = (x[iid]-s.xwin_image)^2+(y[iid]-s.ywin_image)^2
        min=min(dist2, m)
        if min lt 5 then begin
           VisMorphCat[m].morph = morf[iid]
           VisMorphCat[m].ttype = t[iid]
        endif
     endfor
  endif else begin
     readcol, filename, name, ra, dec, ttype, format='A,A,A,I,X'
     for iid=0, n_elements(name)-1 do begin
        rasplit = strsplit(ra[iid], ':', /extract)
        decsplit = strsplit(dec[iid], ':', /extract)
        ra0 = ten(double(rasplit))*15.
        dec0 = ten(double(decsplit))
        w=self->nearestobjs(ra0, dec0)
        VisMorphcat[w[0]].morph = str(ttype[iid])
        VisMorphcat[w[0]].ttype = ttype[iid]
     endfor
  endelse
  w=where(VisMorphCat.ttype ne 99)
  return, VisMorphCat[w]
end

;;;;;;;;;;;;;;;;;;;;;;
;;  Cleanup
Pro Cluster3::CleanUp
  self->FreeImages
  if ptr_valid(self.galapagoscat) then begin
     ptr_free, (*self.galapagoscat).maskfit
     ptr_free, (*self.galapagoscat).simfit
     ptr_free, (*self.galapagoscat).maskfit2
     ptr_free, (*self.galapagoscat).simfit2
  endif
  resolve_routine, 'spectrum__define', /com, /no
  resolve_routine, 'dopplershift__define', /com, /no
  if ptr_valid(self.spectra) then obj_destroy, (*self.spectra).spectrum
  if ptr_valid(self.psf) then ptr_free, *self.psf

  ptr_free, self.band
  ptr_free, self.filehead
  ptr_free, self.zeropoint

  ptr_free, self.SExcat
  ptr_free, self.galapagoscat
  ptr_free, self.skycat
  ptr_free, self.galfitcat
  ptr_free, self.morphcat
  ptr_free, self.morph2cat
  ptr_free, self.CLEANcat
  ptr_free, self.xPSFcat
  ptr_free, self.speccat
  ptr_free, self.spectra

  ptr_free, self.sciim
  ptr_free, self.errim
  ptr_free, self.whtim
  ptr_free, self.skyim
  ptr_free, self.cvlim
  ptr_free, self.mskim
  ptr_free, self.nxpim
  ptr_free, self.hdr
  ptr_free, self.psf
end


;;;;;;;;;;;;;;;;;;;;;;
;;  Define
;;Object definition
Pro Cluster3__Define, struct
  struct = { Cluster3, $
             Inherits JEMobject, $
             rootdir:       '', $
             srcdir:        '', $
             clustername:   '', $
             clusterid:     '', $
             zcluster:      !values.d_nan, $
             veldisp:       !values.d_nan, $
             ra:            !values.d_nan, $
             dec:           !values.d_nan, $
             ebv:           !values.d_nan, $
             platescale:    !values.d_nan, $
             band:          ptr_new(), $
             filehead:      ptr_new(), $
             zeropoint:     ptr_new(), $
             SExcat:        ptr_new(), $
             galapagoscat:  ptr_new(), $
             skycat:        ptr_new(), $
             galfitcat:     ptr_new(), $
             morphcat:      ptr_new(), $
             morph2cat:     ptr_new(), $
             CLEANcat:      ptr_new(), $
             xPSFcat:       ptr_new(), $
             speccat:       ptr_new(), $
             spectra:       ptr_new(), $
             sciim:         ptr_new(), $
             errim:         ptr_new(), $
             whtim:         ptr_new(), $
             skyim:         ptr_new(), $
             cvlim:         ptr_new(), $
             mskim:         ptr_new(), $
             nxpim:         ptr_new(), $
             seghot:        ptr_new(), $
             segcld:        ptr_new(), $
             hdr:           ptr_new(), $
             psf:           ptr_new(), $
             samplefactor:  0L }
end
