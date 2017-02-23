Function Cluster2::Init, $
   _extra=extra

  self->FreeImages
  self->SetProperty, _extra=extra
  return, 1
end

Function Cluster2::SExtract, $
   SExfile, paramfile, band, $
   outdir=outdir, $
   _extra=extra

  if n_elements(outdir) eq 0 then cd, current=outdir
  outdir = directoryify(outdir)
  file_mkdir, outdir
  
  scratchdir = '/scratch-local/jmeyers314/sex'+self.clusterid+'/'
  file_mkdir, scratchdir

  if band eq (*self.band)[0] then begin ;;single image mode
     zeropoint = (*self.zeropoint)[0]
     strsplit = strsplit((*self.filehead)[0], '/', /extract)
     name = strsplit[n_elements(strsplit)-1]

     cmd = 'gunzip -c '+(*self.filehead)[0]+'_SCI.fits.gz > '+scratchdir+name+'_SCI.fits'
     print, cmd
     spawn, cmd
     cmd = 'gunzip -c '+(*self.filehead)[0]+'_ERR.fits.gz > '+scratchdir+name+'_ERR.fits'
     print, cmd
     spawn, cmd

     imagefile = scratchdir+name+'_SCI.fits'
     weightfile = scratchdir+name+'_ERR.fits'
  endif else begin ;;dual image mode
     iband = (where(band eq (*self.band)))[0]
     zeropoint = [(*self.zeropoint)[0], $
                  (*self.zeropoint)[iband]]

     strsplit = strsplit((*self.filehead)[0], '/', /extract)
     name = strsplit[n_elements(strsplit)-1]
     cmd = 'gunzip -c '+(*self.filehead)[0]+'_SCI.fits.gz > '+scratchdir+name+'_SCI.fits'
     print, cmd
     spawn, cmd
     cmd = 'gunzip -c '+(*self.filehead)[0]+'_ERR.fits.gz > '+scratchdir+name+'_ERR.fits'
     print, cmd
     spawn, cmd

     strsplit = strsplit((*self.filehead)[iband], '/', /extract)
     name2 = strsplit[n_elements(strsplit)-1]
     cmd = 'gunzip -c '+(*self.filehead)[iband]+'_SCI.fits.gz > '+scratchdir+name2+'_SCI.fits'
     print, cmd
     spawn, cmd
     cmd = 'gunzip -c '+(*self.filehead)[iband]+'_ERR.fits.gz > '+scratchdir+name2+'_ERR.fits'
     print, cmd
     spawn, cmd
  
     imagefile = [scratchdir+name+'_SCI.fits', $
                  scratchdir+name2+'_SCI.fits']
     weightfile = [scratchdir+name+'_ERR.fits', $
                   scratchdir+name2+'_ERR.fits']
  endelse
  SExtract1, SExfile, paramfile, $
             imagefile, zeropoint, outdir+band+'cat.fits', $
             weightfile=weightfile, weighttype='MAP_RMS', _extra=extra
  cat = mrdfits(outdir+band+'cat.fits',2)
  cmd = 'rm -rf /scratch-local/jmeyers314/sex'+self.clusterid+'/'
  spawn, cmd
  nobjs = (size(cat, /dim))[0]
  galid = lindgen(nobjs)
  s = {clusterid:'', galid:0L, band:band}
  s = replicate(s, nobjs) 
  s.galid=galid
  s.clusterid = self.clusterid
  cat = struct_addtags(s, struct_trimtags(cat,except_tags='NUMBER'))
  return, cat
end

Function Cluster2::SExtractColdHot, $
   coldSExfile, hotSExfile, paramfile, $
   outdir=outdir, $
   store=store
  
  if n_elements(outdir) eq 0 then cd, current=outdir
  outdir=directoryify(outdir)
  ;;do the cold SExtraction
  coldcat = self->SExtract( coldSExfile, paramfile, (*self.band)[0], $
                            outdir=outdir, segment=outdir+'cold_segment.fits' )
  for iband=1, n_elements(*self.band)-1 do begin
     coldcat = [[coldcat], [self->SExtract( coldSExfile, $
                                            paramfile, $
                                            (*self.band)[iband], $
                                            outdir=outdir )]]
  endfor
  nbands = n_elements(*self.band)
  for iband=0, nbands-1 do begin
     filehead = outdir+(*self.band)[iband]+'cat'
     file_copy, filehead+'.fits', filehead+'_cold.fits', /overwrite
  endfor

  ;;now do the hot SExtraction 
  hotcat = self->SExtract( hotSExfile, paramfile, (*self.band)[0], $
                           outdir=outdir, segment=outdir+'hot_segment.fits' )
  for iband=1, n_elements(*self.band)-1 do begin
     hotcat = [[hotcat], [self->SExtract( hotSExfile, $
                                          paramfile, $
                                          (*self.band)[iband], $
                                          outdir=outdir )]]
  endfor
  for iband=0, nbands-1 do begin
     filehead = outdir+(*self.band)[iband]+'cat'
     file_copy, filehead+'.fits', filehead+'_hot.fits', /overwrite
  endfor

  segment = mrdfits(outdir+'cold_segment.fits')
  segment = segment ne 0
  for ix=-2,2 do begin
     segment = segment or shift(segment, ix, 0)
  endfor
  for iy=-2,2 do begin
     segment = segment or shift(segment, 0, iy)
  endfor
  nhotobjs = (size(hotcat, /dim))[0]
  ncoldobjs = (size(coldcat, /dim))[0]
  nbands = (size(coldcat, /dim))[1]
  coldcat = struct_addtags(coldcat, replicate({cold:1}, ncoldobjs, nbands))
  hotcat = struct_addtags(hotcat, replicate({cold:0}, nhotobjs, nbands))
  cat=coldcat
  xsize = (size(segment,/dim))[0]
  ysize = (size(segment,/dim))[1]
  for ihot = 0, nhotobjs-1 do begin
     hotcat1 = hotcat[ihot,*]
     x=fix(hotcat1[0].xwin_image)-1 ;;SEx2IDL coord transform
     y=fix(hotcat1[0].ywin_image)-1
     if x lt 0 or x ge xsize or y lt 0 or y ge ysize then continue
     if segment[x,y] eq 0 then cat = [cat, hotcat1]
  endfor
  nobjs = (size(cat, /dim))[0]
  nbands = (size(cat, /dim))[1]
  galid = lindgen(nobjs)
  cat.galid = rebin(galid, nobjs, nbands)
  if keyword_set(store) then begin
     ptr_free, self.SExcat
     self.SExcat = ptr_new(cat)
  endif
  return, cat
end

Function Cluster2::Galapagos, store=store
  if ~ptr_valid(self.SExcat) then begin
     message, 'Cannot run galapagos without storing a SExtractor catalog', /continue
     return, 0
  endif
  galapagoscat = galapagos1( *self.SExcat )
  if keyword_set(store) then begin
     ptr_free, self.galapagoscat
     self.galapagoscat = ptr_new(galapagoscat)
  endif
  return, galapagoscat
end

Pro Cluster2::LoadImage, band=band
  nbands = n_elements(*self.band)
  if n_elements(band) eq 0 then begin ;;recursively load all images
     for iband=0, nbands-1 do begin
        self->LoadImage, band=(*self.band)[iband]
     endfor
     return
  endif

  if ~ptr_valid(self.image) then begin
     self.image = ptr_new(ptrarr(nbands))
     self.errim = ptr_new(ptrarr(nbands))
     self.hdr = ptr_new(ptrarr(nbands))
  endif

  iband = (where(band eq *self.band))[0]
  ptr_free, (*self.image)[iband], (*self.errim)[iband], (*self.hdr)[iband]

  (*self.image)[iband] = ptr_new(mrdfits((*self.filehead)[iband]+'_SCI.fits.gz', 0))
  (*self.errim)[iband] = ptr_new(mrdfits((*self.filehead)[iband]+'_ERR.fits.gz', 0, hdr1))
  (*self.hdr)[iband] = ptr_new(hdr1)
end

Pro Cluster2::LoadSkyImage, band=band
  nbands = n_elements(*self.band)
  if n_elements(band) eq 0 then begin ;;load all skyimages
     for iband=0, nbands-1 do begin
        self->LoadSkyImage, band=(*self.band)[iband]
     endfor
     return
  endif

  if ~ptr_valid(self.skyim) then begin
     self.skyim = ptr_new(ptrarr(nbands))
  endif

  iband = (where(band eq *self.band))[0]
  ptr_free, (*self.skyim)[iband]

  (*self.skyim)[iband] = ptr_new(mrdfits((*self.filehead)[iband]+'_SKY.fits.gz',0))
end

Pro Cluster2::LoadCVLImage, band=band
  nbands = n_elements(*self.band)
  if n_elements(band) eq 0 then begin ;;load all convolved images
     for iband=0, nbands-1 do begin
        self->LoadCVLImage, band=(*self.band)[iband]
     endfor
     return
  endif

  if ~ptr_valid(self.cvlim) then begin
     self.cvlim = ptr_new(ptrarr(nbands))
  endif

  iband = (where(band eq *self.band))[0]
  ptr_free, (*self.cvlim)[iband]

  (*self.cvlim)[iband] = ptr_new(mrdfits((*self.filehead)[iband]+'_CVL.fits.gz',0))
end

Pro Cluster2::LoadMSKImage, band=band
  nbands = n_elements(*self.band)
  if n_elements(band) eq 0 then begin ;;load all image masks
     for iband=0, nbands-1 do begin
        self->LoadMSKImage, band=(*self.band)[iband]
     endfor
     return
  endif

  if ~ptr_valid(self.mskim) then begin
     self.mskim = ptr_new(ptrarr(nbands))
  endif

  iband = (where(band eq *self.band))[0]
  ptr_free, (*self.mskim)[iband]
  (*self.mskim)[iband] = ptr_new(mrdfits((*self.filehead)[iband]+'_MSK.fits.gz'))
end

Pro Cluster2::LoadWHTImage, band=band
  nbands = n_elements(*self.band)
  if n_elements(band) eq 0 then begin ;;load all image masks
     for iband=0, nbands-1 do begin
        self->LoadWHTImage, band=(*self.band)[iband]
     endfor
     return
  endif

  if ~ptr_valid(self.whtim) then begin
     self.whtim = ptr_new(ptrarr(nbands))
  endif

  iband = (where(band eq *self.band))[0]
  ptr_free, (*self.whtim)[iband]
  (*self.whtim)[iband] = ptr_new(mrdfits((*self.filehead)[iband]+'_WHT.fits.gz'))
end

Pro Cluster2::LoadSegHotImage
  id=self.clusterid
  ptr_free, self.seghot
  self.seghot = ptr_new(mrdfits('/home/scpdata02/clusters/'+id+'/hot_segment.fits'))
end

Pro Cluster2::LoadSegColdImage
  id=self.clusterid
  ptr_free, self.segcold
  self.segcold = ptr_new(mrdfits('/home/scpdata02/clusters/'+id+'/cold_segment.fits'))
end

Pro Cluster2::FreeImages
  if ptr_valid(self.image) then ptr_free, *self.image
  if ptr_valid(self.errim) then ptr_free, *self.errim
  if ptr_valid(self.skyim) then ptr_free, *self.skyim
  if ptr_valid(self.cvlim) then ptr_free, *self.cvlim
  if ptr_valid(self.mskim) then ptr_free, *self.mskim
  if ptr_valid(self.whtim) then ptr_free, *self.whtim
  ptr_free, self.seghot
  ptr_free, self.segcold
  if ptr_valid(self.hdr) then ptr_free, *self.hdr
end

Function Cluster2::ImageSection, $
   p, $      ;; position to extract
   band=band, $
   err=err, $
   cvlim=cvlim, $
   segcold=segcold, $
   seghot=seghot

  if n_elements(band) eq 0 then band=(*self.band)[0]
  iband = (where(band eq *self.band))[0]
  if ~ptr_valid((*self.image)[iband]) then begin
     out = dblarr(p[2]-p[0]+1,p[3]-p[1]+1)*!values.d_nan
  endif else begin
     sz=size((*(*self.image)[iband]),/dim)
     if p[0] lt 0. or p[2] lt 0. or p[1] gt sz[0]-1 or p[3] gt sz[1]-1 then return, $
        dblarr(p[2]-p[0]+1,p[3]-p[1]+1)*!values.d_nan
     out = (*(*self.image)[iband])[p[0]:p[2],p[1]:p[3]]
  endelse
  if arg_present(err) then $
     err = (*(*self.errim)[iband])[p[0]:p[2],p[1]:p[3]]
  if arg_present(cvlim) then $
     cvlim = (*(*self.cvlim)[iband])[p[0]:p[2],p[1]:p[3]]
  if arg_present(segcold) then $
     segcold = (*self.segcold)[ p[0]:p[2], p[1]:p[3] ]
  if arg_present(seghot) then $
     seghot = (*self.seghot)[ p[0]:p[2], p[1]:p[3] ]
  return, out
end

Function Cluster2::GalImage, gal, $
                             xsize=xsize, $
                             ysize=ysize, $
                             skysub=skysub, $
                             err=err, $
                             _extra=extra
  if n_elements(xsize) eq 0 or n_elements(ysize) eq 0 then begin
     gal1=(*self.galapagoscat)[gal]
     p=[gal1.xmin,gal1.ymin,gal1.xmax,gal1.ymax]
  endif else begin
     s=(*self.SExcat)[gal,0]
     x=s.xwin_image
     y=s.ywin_image
     p=[x-xsize/2, y-ysize/2, x+xsize/2, y+ysize/2]
  endelse
  if keyword_set(err) then junk = self->ImageSection(p, err=img, _extra=extra) $
  else img = self->ImageSection(p, _extra=extra)
  if keyword_set(skysub) then begin
     if n_elements(band) eq 0 then band = (*self.band)[0]
     skyvalue = (*self.skycat)[gal,where(band eq *self.band)].sky_value
     if finite(skyvalue) then img -= skyvalue
  endif
  return, img
end

Function Cluster2::NewIsophotSky, gals=gals, band=band, store=store
  if ~ptr_valid(self.SExcat) || ~ptr_valid(self.galapagoscat) || $
     ~ptr_valid(self.errim) || ~ptr_valid(self.image) || ~ptr_valid(self.whtim) then begin
     message, 'IsophotSky requires a stored SExtractor catalog and stored GALAPAGOS catalog', /continue
     return, 0
  endif
  if n_elements(gals) eq 0 then gals=indgen(n_elements(galapagoscat))
  store=keyword_set(store)

  if ~ptr_valid(self.skycat) then begin
     skycat = { skycat, $
                clusterid: '', $
                galid:     0L, $
                band:      '', $
                sky_value: !values.d_nan, $
                sky_sigma: !values.d_nan, $
                sky_total: !values.d_nan, $
                exptime: !values.d_nan }
     nobjs = (size(*self.SExcat, /dim))[0]
     nbands = (size(*self.SExcat, /dim))[1]
     skycat = replicate( skycat, nobjs, nbands )
     skycat.clusterid = self.clusterid
     skycat.galid = rebin( lindgen(nobjs), nobjs, nbands )
     for iband=0, n_elements(*self.band)-1 do begin
        skycat[*,iband].band=(*self.band)[iband]
     endfor
     self.skycat = ptr_new(skycat)
  endif

  skycat = (*self.skycat)  
  if n_elements(band) eq 0 then begin
     for iband=0, n_elements(*self.band)-1 do begin
        skycat[gals,iband] = (self->NewIsophotSky(gals=gals, $
                                                  band=(*self.band)[iband], $
                                                  store=store ))[*,iband]
     endfor
     return, skycat
  endif

  iband = (where(band eq (*self.band)))[0]
  platescale = sxpar(*(*self.hdr)[iband], 'D001SCAL')*0.05
  galapagoscat = *self.galapagoscat
  SExcat = *self.SExcat
  @definestructs
  platefrac = 0.05/platescale
  sz = size(*(*self.image)[iband], /dim)
  for igal=0, n_elements(gals)-1 do begin
     counter, igal+1, n_elements(gals), 'IsophotSky on band '+band+': '
     gal1=galapagoscat[gals[igal]]
     sex1=SExcat[gals[igal],0]
     ;; check if center is off the edge of all exposures...
    
     xmin = round(gal1.xmin - 250*platefrac > 0)
     xmax = round(gal1.xmax + 250*platefrac < (sz[0]-1))
     ymin = round(gal1.ymin - 250*platefrac > 0)
     ymax = round(gal1.ymax + 250*platefrac < (sz[1]-1))
     sciim = (*(*self.image)[iband])[xmin:xmax, ymin:ymax]
     errim = (*(*self.errim)[iband])[xmin:xmax, ymin:ymax]
     mskim = (*(*self.mskim)[iband])[xmin:xmax, ymin:ymax]
     ellipse = ellipse0
     ellipse.center = [sex1.xwin_image-xmin-1, sex1.ywin_image-ymin-1]
     ellipse.majaxis = sex1.a_image*sex1.kron_radius
     ellipse.minaxis = sex1.b_image*sex1.kron_radius
     ellipse.theta = sex1.theta_image
     skystruct = newisophotsky( sciim, mskim, ellipse, platescale=platescale )
     skycat[gals[igal], iband].sky_value = skystruct.sky_value
     skycat[gals[igal], iband].sky_sigma = skystruct.sky_sigma
     majaxis=sex1.a_image*sex1.kron_radius
     xmin = gal1.xmin > 0
     ymin = gal1.ymin > 0
     xmax = gal1.xmax < (sz[0]-1)
     ymax = gal1.ymax < (sz[1]-1)
     skyim = (*(*self.skyim)[iband])[xmin:xmax, ymin:ymax]
     whtim = (*(*self.whtim)[iband])[xmin:xmax, ymin:ymax]
     dist_ellipse, mask, [xmax-xmin+1, ymax-ymin+1], sex1.xwin_image-xmin-1, sex1.ywin_image-ymin-1, sex1.a_image/sex1.b_image, sex1.theta_image-90.
     e=mask le majaxis
     if total(e) eq 0 then begin
        skycat[gals[igal], iband].sky_total = skyim[(xmax-xmin+1)/2., (ymax-ymin+1)/2.]
        skycat[gals[igal], iband].exptime = whtim[(xmax-xmin+1)/2., (ymax-ymin+1)/2.]
     endif else begin
        skycat[gals[igal], iband].sky_total = biweight_mean((skyim*e)[where(skyim*e ne 0.)])
        skycat[gals[igal], iband].exptime = biweight_mean((whtim*e)[where(skyim*e ne 0.)])
     endelse
  endfor
  print
  if keyword_set(store) then begin
     ptr_free, self.skycat
     self.skycat = ptr_new(skycat)
  endif
  return, skycat[gals,*]
end

Function Cluster2::IsophotSky, gals=gals, band=band, store=store
  if ~ptr_valid(self.SExcat) || ~ptr_valid(self.galapagoscat) then begin
     message, 'IsophotSky requires a stored SExtractor catalog and stored GALAPAGOS catalog', /continue
     return, 0
  endif
  if n_elements(gals) eq 0 then gals=indgen(n_elements(galapagoscat))
  store=keyword_set(store)

  if ~ptr_valid(self.skycat) then begin
     skycat = { skycat, $
                clusterid: '', $
                galid:     0L, $
                band:      '', $
                sky_value: 0.d, $
                sky_sigma: 0.d, $
                sky_total: 0.d }
     nobjs = (size(*self.SExcat, /dim))[0]
     nbands = (size(*self.SExcat, /dim))[1]
     skycat = replicate( skycat, nobjs, nbands )
     skycat.clusterid = self.clusterid
     skycat.galid = rebin( lindgen(nobjs), nobjs, nbands )
     for iband=0, n_elements(*self.band)-1 do begin
        skycat[*,iband].band=(*self.band)[iband]
     endfor
     self.skycat = ptr_new(skycat)
  endif

  skycat = (*self.skycat)  
  if n_elements(band) eq 0 then begin
     for iband=0, n_elements(*self.band)-1 do begin
        skycat[gals,iband] = (self->IsophotSky(gals=gals, $
                                               band=(*self.band)[iband], $
                                               store=store ))[*,iband]
     endfor
     return, skycat
  endif

  iband = (where(band eq (*self.band)))[0]
  platescale = sxpar(*(*self.hdr)[iband], 'D001SCAL')*0.05
  galapagoscat = *self.galapagoscat
  SExcat = *self.SExcat

  for igal=0, n_elements(gals)-1 do begin
     counter, igal+1, n_elements(gals), 'IsophotSky on band '+band+': '
     IsophotSky, (*self.image)[iband], (*self.errim)[iband], $
                 SExcat, galapagoscat, $
                 gals[igal], $
                 platescale=platescale, $
                 sky_value=sky_value, $
                 sky_sigma=sky_sigma
     skycat[gals[igal],iband].sky_value = sky_value
     skycat[gals[igal],iband].sky_sigma = sky_sigma
     s0 = SExcat[gals[igal],0]
     x = s0.xwin_image
     y = s0.ywin_image
     majaxis = s0.a_image*s0.kron_radius
     xmin = floor(x-majaxis-10)
     ymin = floor(y-majaxis-10)
     xmax = floor(x+majaxis+10)
     ymax = floor(y+majaxis+10)
     sz=size((*(*self.skyim)[iband]), /dim)
     xmin = xmin > 0
     ymin = ymin > 0
     xmax = xmax < (sz[0]-1)
     ymax = ymax < (sz[1]-1)
     if xmin gt xmax or ymin gt ymax then begin
        skycat[gals[igal],iband].sky_total = !values.f_nan
     endif else begin
        skyim = (*(*self.skyim)[iband])[xmin:xmax, ymin:ymax]
        dist_ellipse, mask, [xmax-xmin+1, ymax-ymin+1], x-xmin-1, y-ymin-1, s0.a_image/s0.b_image, s0.theta_image-90.
        e=mask le majaxis
        if total(e) eq 0 then skycat[gals[igal],iband].sky_total = skyim[xmax-xmin+1, ymax-ymin+1] $
        else skycat[gals[igal],iband].sky_total = biweight_mean((skyim*e)[where(skyim*e ne 0.)])
     endelse
  endfor
  print
  if keyword_set(store) then begin
     ptr_free, self.skycat
     self.skycat = ptr_new(skycat)
  endif
  return, skycat[gals,*]
end

Function Cluster2::IsophotskyPlot, gal, band=band, $
                                   annuli=annuli, mask=mask, sky=sky, grad=grad, finalmask=finalmask
  if n_elements(band) eq 0 then band = (*self.band)[0]
  iband = (where(band eq *self.band))[0]

  gal1=(*self.galapagoscat)[gal]
  SEx1=(*self.SExcat)[gal,iband]
  
  platescale = sxpar(*(*self.hdr)[0], 'D001SCAL')*0.05
  platefrac = 0.05/platescale
  xmin = round(gal1.xmin - 250*platefrac > 0)
  xmax = round(gal1.xmax + 250*platefrac < (size((*(*self.image)[iband]), /dim))[0]-1)
  ymin = round(gal1.ymin - 250*platefrac > 0)
  ymax = round(gal1.ymax + 250*platefrac < (size((*(*self.image)[iband]), /dim))[1]-1)

  nx = xmax - xmin + 1
  ny = ymax - ymin + 1
  
  skyim = (*(*self.image)[iband])[xmin:xmax, ymin:ymax]
  skyer = (*(*self.errim)[iband])[xmin:xmax, ymin:ymax]
  mask  = skyim*0

  if ~finite((*(*self.errim)[iband])[(gal1.xmin+gal1.xmax)/2,(gal1.ymin+gal1.ymax)/2]) $
  then return, 0
  
  ellipses=bytarr( nx, ny, 14 )
  annuli  =bytarr( nx, ny, 12 )
  dist_ellipse, skyellipse, [nx,ny], $
                SEx1.xwin_image-xmin-1, SEx1.ywin_image-ymin-1, $
                SEx1.a_image/SEx1.b_image, SEx1.theta_image-90.
  
  for iellipse=0, 13 do begin
     ellipses[*,*,iellipse] = skyellipse lt gal1.majaxis + 18*iellipse*platefrac
  endfor

  for iannulus=0, 11 do begin
     annuli[*,*,iannulus] = ellipses[*,*,iannulus+2] and (1 - ellipses[*,*,iannulus])
  endfor

  w = where( (*self.galapagoscat).xmin le xmax and $
             (*self.galapagoscat).xmax ge xmin and $
             (*self.galapagoscat).ymin le ymax and $
             (*self.galapagoscat).ymax ge ymin )
  if w[0] ne -1 then begin
     for j=0, n_elements(w)-1 do begin
        maskgal = (*self.galapagoscat)[w[j]]
        maskSEx = (*self.SExcat)[w[j],0]
        dist_ellipse, maskellipse, [nx,ny], $
                      maskSEx.xwin_image-xmin-1, maskSEx.ywin_image-ymin-1, $
                      maskSEx.a_image/maskSEx.b_image, maskSEx.theta_image-90.
        e = maskellipse le maskgal.majaxis
        mask = mask or e
     endfor
  endif
  mask = mask or (1-finite( skyer ))
  
  sky = dblarr(12)
  minstart=0
  for iannulus=0, 11 do begin
     wsky = where( annuli[*,*,iannulus] and (1-mask) )
     if n_elements(wsky) lt 20 then minstart += 1 $
     else sky[iannulus] = biweight_mean( skyim[wsky] )
  endfor

  grad=dblarr(7)
  for igrad=0,6 do begin
     result = linfit(indgen(6), sky[igrad:igrad+5])
     grad[igrad]=result[1]
  endfor

  wgrad = where(grad gt -0.003/(18*platefrac)*0.01)
  if wgrad[0] eq -1 then istart=6 $
  else istart=max([min(wgrad),minstart]) < 6
  finalmask=ellipses[*,*,istart+6] and (1-ellipses[*,*,istart]) and (1-mask)

  return, skyim
end

Function Cluster2::GetStars, band=band, starsize=starsize, $
                             rthresh=rthresh, minstars=minstars
  if n_elements(band) eq 0 then band=(*self.band)[0]
  iband = (where(band eq *self.band))[0]
  SExcat = *self.SExcat
  skycat = *self.skycat
  platefrac = sxpar(*(*self.hdr)[0], 'D001SCAL')
  platescale = platefrac*0.05
  if n_elements(starsize) eq 0 then begin
     starsize=fix(41/platefrac)
     if fix(starsize/2.) eq starsize/2. then starsize+=1 ;;ensure odd
  endif
  innersize = fix(starsize*0.23)
  if fix(innersize/2.) eq innersize/2. then innersize+=1 ;;ensure odd

  w = where( SExcat[*,0].flux_radius gt 1.4/platefrac and $
             SExcat[*,0].flux_radius lt 1.72/platefrac and $
             SExcat[*,0].mag_auto gt 18 and $
             SExcat[*,0].mag_auto lt 24 ) 

  if n_elements(rthresh) ne 0 then begin
     ;; radial threshold in arcmin.
     ra = SExcat[w,0].alphawin_j2000
     dec = SExcat[w,0].deltawin_j2000
     ra0 = self.ra
     dec0 = self.dec
     r2 = (ra0-ra)^2*(cos(dec0*!dpi/180.d))^2 $
          + (dec0-dec)^2
     r = sqrt(r2)*60. ;convert degree->arcmin
     if n_elements(minstars) eq 0 then begin
        w=w[where(r le rthresh)]
     endif else begin
        if n_elements(w) lt minstars then begin
           print, 'not enough stars!'
           print, 'using all stars found!'
        endif else begin
           if total(r le rthresh) lt minstars then begin
              s=sort(r)
              w=w[s[0:minstars-1]]
           endif else begin
              w=w[where(r le rthresh)]
           endelse
        endelse
     endelse
  endif
     

  halfsize=(starsize-1)/2
  onestar = {clusterid:'', $
             galid:0L, $
             image:dblarr(starsize,starsize), $
             model:dblarr(innersize,innersize), $
             errim:dblarr(starsize,starsize), $
             modelparams: dblarr(7), $
             x_SEx:0.d, $
             y_SEx:0.d }
  for iw=0, n_elements(w)-1 do begin
     star = onestar
     star.clusterid=self.clusterid
     star.galid=w[iw]
     p = round([SExcat[w[iw],0].xwin_image-1, $ ;;need to take SEx coords -> IDL coords
                SExcat[w[iw],0].ywin_image-1, $
                SExcat[w[iw],0].xwin_image-1, $
                SExcat[w[iw],0].ywin_image-1])+[-1,-1,1,1]*halfsize
     star.image = self->ImageSection(p, band=band, err=err)
     if n_elements(err) eq 0 || total(~finite(err)) ne 0 then continue
     sky_value = skycat[w[iw],iband].sky_value
     if finite(sky_value) then star.image -= sky_value
     star.errim = err
     star.x_SEx = SExcat[w[iw],0].xwin_image-fix(SExcat[w[iw],0].xwin_image)
     star.y_SEx = SExcat[w[iw],0].ywin_image-fix(SExcat[w[iw],0].ywin_image)

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

Function Cluster2::Galfit, $
   gals=gals, $
   outdir=outdir, $
   take2=take2, $
   noerr=noerr, $
   store=store, $
   doall=doall

  if ~ptr_valid(self.galapagoscat) || ~ptr_valid(self.SExcat) $
     || ~ptr_valid(self.skycat) || ~ptr_valid(self.image) $
     || ~ptr_valid((*self.image)[0]) || ~ptr_valid(self.psf) $
     || ~ptr_valid((*self.psf)[0]) then begin
     message, 'Galfit requires a stored SExtractor catalog, a stored GALAPAGOS catalog,' $
              + ' a stored sky background catalog, and stored images/PSFs', /continue
     return, 0
  endif

  if n_elements(outdir) eq 0 then cd, current=outdir
  if n_elements(gals) eq 0 then gals = lindgen( n_elements( *self.galapagoscat ) )
  if total(gals) eq -1 then return, 0
  take2 = keyword_set(take2)
  noerr = keyword_set(noerr)
  doall = keyword_set(doall)
  store = keyword_set(store)

  if ~ptr_valid(self.galfitcat) then begin
     galfitcat1 = { galfitcat, $
                    clusterid:  '', $
                    galid:      0L, $
                    xpos:       -1.d, $
                    ypos:       -1.d, $
                    mag:        -1.e, $
                    Re:         -1.e, $
                    boa:        -1.e, $
                    pa:         -1.e, $
                    n:          -1.e, $
                    errxpos:    -1.e, $
                    errypos:    -1.e, $
                    errmag:     -1.e, $
                    errRe:      -1.e, $
                    errboa:     -1.e, $
                    errpa:      -1.e, $
                    errn:       -1.e, $
                    chi2:       -1.e, $
                    dof:        -1, $
                    chi2perdof: -1.e, $
                    galfittime: -1.e, $
                    galfitretry: -1 }
     galfitcat = replicate( galfitcat1, n_elements(*self.galapagoscat) )
     galfitcat.clusterid = (*self.galapagoscat).clusterid
     galfitcat.galid = (*self.galapagoscat).galid
     self.galfitcat = ptr_new(galfitcat)
  endif

  if doall then begin
     galfitcat0 = self->Galfit( gals=gals, outdir=outdir, store=store )
     s=self->summary()
     wfail = where( s.galfittime ne -1 and s.re eq -1 )
     galfitcat1 = self->Galfit( gals=wfail, outdir=outdir, /take2, store=store )
     s=self->summary()
     wfail2 = where( s.galfittime ne -1 and s.re eq -1 )
     galfitcat2 = self->Galfit( gals=wfail2, outdir=outdir, /noerr, store=store )
     s=self->summary()
     wfail3 = where( s.galfittime ne -1 and s.re eq -1 )
     galfitcat3 = self->Galfit( gals=wfail3, outdir=outdir, /take2, /noerr, store=store )
     if wfail[0] ne -1 then galfitcat0[wfail]=galfitcat1
     if wfail2[0] ne -1 then galfitcat0[wfail2]=galfitcat2
     if wfail3[0] ne -1 then galfitcat0[wfail3]=galfitcat3
     return, galfitcat0
  endif else begin
     
     galfitcat = (*self.galfitcat)
     psf = undersampleimage(*(*self.psf)[0], self.samplefactor)
     psf /= max(psf)
     psf_ptr = ptr_new(psf)
     for igal=0L, n_elements(gals)-1 do begin
        starttime = systime(1)
        galap1 = (*self.galapagoscat)[gals[igal]]
        mkhdr, hdr, 4, [1,1]
        sxaddpar, hdr, 'EXPTIME', (*self.skycat)[gals[igal],0].exptime
        sxaddpar, hdr, 'GAIN', 1.0
        sxaddpar, hdr, 'RDNOISE', 5.0
        sxaddpar, hdr, 'NCOMBINE', sxpar(*(*self.hdr)[0], 'NDRIZIM')/2
        galfit1 = MakeGalfit1( (*self.image)[0], (*self.errim)[0], psf_ptr, $
                               *self.galapagoscat, *self.SExcat, *self.skycat, $
                               gals[igal], (*self.band)[0], 0, $
                               (*self.zeropoint)[0], (*self.skycat)[gals[igal],0].sky_total, $
                               outdir=outdir, header=hdr, take2=take2, noerr=noerr)
        if size(galfit1, /tname) ne 'OBJREF' then continue
        galfit1->Execute
        results = galfit1->Results()
        if size(results, /tname) eq 'STRUCT' then begin
           results.xpos += galap1.xmin
           results.ypos += galap1.ymin
           dest = galfitcat[gals[igal]]
           struct_assign, results, dest, /nozero
           galfitcat[gals[igal]] = dest
           case 1 of 
              (~keyword_set(take2) and ~keyword_set(noerr)) : galfitcat[gals[igal]].galfitretry=0
              (keyword_set(take2) and ~keyword_set(noerr))  : galfitcat[gals[igal]].galfitretry=1
              (~keyword_set(take2) and keyword_set(noerr))  : galfitcat[gals[igal]].galfitretry=2
              (keyword_set(take2) and keyword_set(noerr))   : galfitcat[gals[igal]].galfitretry=3
           endcase        
        endif
        galfitcat[gals[igal]].galfittime = systime(1)-starttime
        obj_destroy, galfit1
     endfor
     ptr_free, psf_ptr
     if keyword_set(store) then begin
        ptr_free, self.galfitcat
        self.galfitcat = ptr_new(galfitcat)
     endif
  endelse
  return, galfitcat[gals]
end

Function Cluster2::GalfitObj, $
   gal, $
   band=band, $
   outdir=outdir, $
   _extra=extra
;; extra includes:
;; noerr=noerr
;; take2=take2 for 2nd chance GALFITing.
;; snx=snx
;; sny=sny
;; snmag=snmag

  if n_elements(band) eq 0 then band=(*self.band)[0]
  iband = (where(band eq (*self.band)))[0]
  if n_elements(outdir) eq 0 then cd, current=outdir
  psf = undersampleimage(*(*self.psf)[iband], self.samplefactor)
  psf /= max(psf)
  psf_ptr = ptr_new(psf)
  header = *(*self.hdr)[iband]
  imageptr = (*self.image)[iband]
  errimptr = (*self.errim)[iband]

  mkhdr, hdr, 4, [1,1]
  sxaddpar, hdr, 'EXPTIME', (*self.skycat)[gal,iband].exptime
  sxaddpar, hdr, 'GAIN', 1.0
  sxaddpar, hdr, 'RDNOISE', 5.0
  sxaddpar, hdr, 'NCOMBINE', sxpar(header, 'NDRIZIM')/2
  galfit1 = MakeGalfit1( imageptr, errimptr, psf_ptr, $
                         *self.galapagoscat, *self.SExcat, *self.skycat, $
                         gal, band, iband, $
                         (*self.zeropoint)[iband], (*self.skycat)[gal,iband].sky_total, $
                         outdir=outdir, header=hdr, _extra=extra )
  ptr_free, psf_ptr
  return, galfit1
end

Function Cluster2::RecallGalfit, gal, outdir=outdir
  if n_elements(outdir) eq 0 then outdir='/home/scpdata02/clusters/'+self.clusterid+'/galfit/'
  outdir = directoryify(outdir)
  prefix=self.clusterid+string(gal,format='(I05)')+'z'
  filename = outdir+prefix+'.fits'
  if ~(file_info(filename)).exists then return, -1.
  struct = mrdfits(filename,1,/silent)
  return, struct
end

Function Cluster2::Morphology, $
   gals=gals, $
   outdir=outdir, $
   store=store
   
  if ~ptr_valid(self.morphcat) then begin
     morphcat = { morphcat, $
                  clusterid: '', $
                  galid:     0L, $
                  qpflux:    !values.f_nan, $
                  qprad:     !values.f_nan, $
                  prad:      !values.f_nan, $
                  asym:      !values.f_nan, $
                  asymerr:   !values.f_nan, $
                  gini:      !values.f_nan, $
                  ginierr:   !values.f_nan, $
                  m20:       !values.f_nan, $
                  conc:      !values.f_nan, $
                  concrank:  !values.f_nan, $
                  smooth:    !values.f_nan, $
                  avsnr:     !values.f_nan }
     morphcat = replicate(morphcat, n_elements(*self.galfitcat))
     morphcat.clusterid = (*self.galfitcat).clusterid
     morphcat.galid = (*self.galfitcat).galid
     self.morphcat = ptr_new(morphcat)
  endif
  morphcat = *self.morphcat
  
  if n_elements(gals) eq 0 then gals=lindgen(n_elements(*self.galapagoscat))
  for igal=0, n_elements(gals)-1 do begin
     counter, igal+1, n_elements(gals), 'Measuring Morphology:'

     a=self->recallgalfit( gals[igal], outdir=outdir )
     if size(a,/tname) ne 'STRUCT' then continue

     skysig = (*self.skycat)[gals[igal],0].sky_sigma
     gal1=(*self.galfitcat)[gals[igal]]
     if member(gal1.galfitretry,[2,3]) then begin
        exptime = (*self.skycat)[gals[igal],0].exptime
        a.intsb /= exptime ;; recall intsb is already sky-subtracted
        a.errim /= exptime
     endif
     center = [a.xpos, a.ypos]-1.  ;galfit2IDL coord shift
     acenter = center
     s=size(a.intsb,/dim)
     if center[0] lt 0. $
        or center[1] lt 0. $
        or center[0] gt (s[0]-1) $
        or center[1] gt (s[1]-1) then continue
     if total(finite(a.intsb)) ne s[0]*s[1] then continue
     
     seg = ACSsegmask( a.intsb, 1.5*skysig )
     mask = seg eq seg[center[0], center[1]]
     if total(mask) eq 0 then continue
     mask = mask or shift(mask, 1)
     mask = mask or shift(mask, -1)
     mask = mask or shift(mask, 0, 1)
     mask = mask or shift(mask, 0, -1)
     prad = petrorad( a.intsb, center, a.boa, a.pa )
     qprad = quasipetrorad( a.intsb, mask, center, a.boa, a.pa, qpflux=qpflux )
     
     if ~finite(qpflux) then continue

     mask = a.intsb*mask gt qpflux
     if total(mask) eq 0 then continue
     avsnr = mean((a.intsb/a.errim)[where(mask)])
     gini = gini( a.intsb[where(mask)], err=ginierr )
     asym = asymmetry( a.intsb, $
                       mask, $
                       skysig, $
                       center=acenter, $
                       err=asymerr, $
                       /interp )
     guess=center
     m20 = m20( a.intsb, mask, guess=guess )
     conc = concentration( a.intsb, center, qprad )
     concrank = concrank( a.intsb[where(mask)] )
     smooth=smoothness( a.intsb, skysig, qprad, center=center )

     morphcat[gals[igal]].qpflux=qpflux
     morphcat[gals[igal]].qprad=qprad
     morphcat[gals[igal]].prad=prad
     morphcat[gals[igal]].gini=gini
     morphcat[gals[igal]].ginierr=ginierr
     morphcat[gals[igal]].asym=asym
     morphcat[gals[igal]].asymerr=asymerr
     morphcat[gals[igal]].m20=m20
     morphcat[gals[igal]].conc=conc
     morphcat[gals[igal]].concrank=concrank
     morphcat[gals[igal]].smooth=smooth
     morphcat[gals[igal]].avsnr=avsnr
  endfor
  if keyword_set(store) then begin
     ptr_free, self.morphcat
     self.morphcat = ptr_new(morphcat)
  endif
  return, morphcat[gals]
end

Function Cluster2::Morph2, $
   gals=gals, $
   outdir=outdir, $
   store=store

  if ~ptr_valid(self.morph2cat) then begin
     morph2cat = { morph2cat, $
                   clusterid: '', $
                   galid:     0L, $
                   qpflux2:    !values.f_nan, $
                   qprad2:     !values.f_nan, $
                   prad2:      !values.f_nan, $
                   asym2:      !values.f_nan, $
                   asymerr2:   !values.f_nan, $
                   gini2:      !values.f_nan, $
                   ginierr2:   !values.f_nan, $
                   m202:       !values.f_nan, $
                   conc2:      !values.f_nan, $
                   concrank2:  !values.f_nan, $
                   smooth2:    !values.f_nan, $
                   avsnr2:     !values.f_nan }
     morph2cat = replicate(morph2cat, n_elements(*self.galfitcat))
     morph2cat.clusterid = (*self.galfitcat).clusterid
     morph2cat.galid = (*self.galfitcat).galid
     self.morph2cat = ptr_new(morph2cat)
  endif
  morph2cat = *self.morph2cat
  
  if n_elements(gals) eq 0 then gals=lindgen(n_elements(*self.galapagoscat))
  for igal=0, n_elements(gals)-1 do begin
     counter, igal+1, n_elements(gals), 'Measuring Morphology:'

     SExcat1 = (*self.SExcat)[gals[igal],0]
     Galcat1 = (*self.galfitcat)[gals[igal]]
     a=self->recallgalfit(gals[igal])
     if size(a, /tname) ne 'STRUCT' then continue
     sz = size(a.image, /dim)
     p=intarr(4)
     p[0] = Galcat1.xpos-a.xpos
     p[1] = Galcat1.ypos-a.ypos
     p[2] = p[0]+sz[0]-1
     p[3] = p[1]+sz[1]-1
     if SExcat1.cold then $
        junk = self->ImageSection( p, band='z', segcold=seg ) $
     else $
        junk = self->ImageSection( p, band='z', seghot=seg )
     a.image -= (*self.skycat)[gals[igal],0].sky_value
     skysig=(*self.skycat)[gals[igal],0].sky_sigma
     center = floor([a.xpos, a.ypos]-1.)
     if center[0] lt 0 or center[0] ge sz[0] or $
        center[1] lt 0 or center[1] ge sz[1] then continue
     acenter=center
     mask = seg eq seg[center[0], center[1]]
     if total(mask) eq 0 then continue
     mask = mask or shift(mask, 1)
     mask = mask or shift(mask, -1)
     mask = mask or shift(mask, 0, 1)
     mask = mask or shift(mask, 0, -1)
     prad = petrorad( a.image, center, a.boa, a.pa )
     qprad = quasipetrorad( a.image, mask, center, a.boa, a.pa, qpflux=qpflux )
     if ~finite(qpflux) then continue
     mask = a.image*mask gt qpflux
     if total(mask) eq 0 then continue
     avsnr = mean((a.image/a.errim)[where(mask)])
     gini = gini( a.image[where(mask)], err=ginierr )
     asym = asymmetry( a.image, $
                       mask, $
                       skysig, $
                       center=acenter, $
                       err=asymerr, $
                       /interp )
     guess=center
     m20 = m20( a.image, mask, guess=guess )
     conc = concentration( a.image, center, qprad )
     concrank = concrank( a.image[where(mask)] )
     smooth=smoothness( a.image, skysig, qprad, center=center )

     morph2cat[gals[igal]].qpflux2=qpflux
     morph2cat[gals[igal]].qprad2=qprad
     morph2cat[gals[igal]].prad2=prad
     morph2cat[gals[igal]].gini2=gini
     morph2cat[gals[igal]].ginierr2=ginierr
     morph2cat[gals[igal]].asym2=asym
     morph2cat[gals[igal]].asymerr2=asymerr
     morph2cat[gals[igal]].m202=m20
     morph2cat[gals[igal]].conc2=conc
     morph2cat[gals[igal]].concrank2=concrank
     morph2cat[gals[igal]].smooth2=smooth
     morph2cat[gals[igal]].avsnr2=avsnr
  endfor
  if keyword_set(store) then begin
     ptr_free, self.morph2cat
     self.morph2cat = ptr_new(morph2cat)
  endif
  return, morph2cat[gals]

end

Function Cluster2::Morph2Plot, gal, outdir=outdir
  acskern = [[0.0322, 0.0718, 0.0322], $
             [0.0718, 0.1600, 0.0718], $
             [0.0322, 0.0718, 0.0322]]
  acskern /= total(acskern)
  skycat1 = (*self.skycat)[gal,0]
  a=self->recallgalfit(gal, outdir=outdir)
  if size(a,/tname) ne 'STRUCT' then return, !values.f_nan
  s=self->summary()
  if member( s[gal].galfitretry, [2,3] ) then begin
     exptime = (*self.skycat)[gal,0].exptime
     a.intsb /= exptime
     a.errim /= exptime
  endif
  s=size(a.intsb,/dim)
  seg = ACSsegmask( a.intsb, 1.5*skycat1.sky_sigma )
  center = [a.xpos, a.ypos]-1.
  if center[0] lt 0. $
     or center[1] lt 0. $
     or center[0] gt (s[0]-1) $
     or center[1] gt (s[1]-1) then return, !values.f_nan
  mask = seg eq seg[center[0], center[1]]
  if total(mask) eq 0 then return, !values.f_nan
  wmask = where( mask )
  qpflux = quasipetroflux( a.intsb*mask )
  w=where(a.intsb*mask gt qpflux)
  if w[0] eq -1 then return, !values.f_nan
  gini = gini( a.intsb[w] )
  asym = asymmetry( a.intsb, $
                    a.intsb*mask gt qpflux, $
                    skycat1.sky_sigma, $
                    center=center, $
                    rotresid=rotresid)
  return, [[mask*a.intsb, a.intsb*(a.intsb*mask gt qpflux), rotresid], $
           [a.image, a.model, a.resid]]
end

Pro Cluster2::AddSpectroscopy, $
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

Pro Cluster2::MeasureO2, plot=plot1
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
              plot='/home/scpdata02/clusters/analysis/o2/'+name+'.o2.eps'
              O2_EW=measureo2(s1, plot=plot, err=O2_EWerr)
           endif else O2_EW=measureo2(s1, err=O2_EWerr)
           spectra[w[iw]].o2_ew=O2_EW
           spectra[w[iw]].o2_ewerr=O2_EWerr
        endif
     endfor
  endfor
  ptr_free, self.spectra
  self.spectra=ptr_new(spectra)
  print
end

Pro Cluster2::UpdateSpecCat
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
  if ptr_valid(self.spectra) then begin
     gals = ((*self.spectra).galid)[uniq( (*self.spectra).galid, $
                                          sort( (*self.spectra).galid ) )]
  endif
  for igal=0, n_elements(gals)-1 do begin
     if gals[igal] eq -1 then continue;; these are the unidentified spectra...
     w = where( (*self.spectra).galid eq gals[igal] )
     wgood = where( (*self.spectra)[w].z ne -1. )
     if wgood[0] eq -1 then begin  ;;didn't get any redshifts of this object
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
end

Pro Cluster2::CreateRegion, filename, gals=gals
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
     gcat = (*self.galapagoscat)[gals[i]]
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

Pro Cluster2::SpecRegionCheck, filename
  if ~ptr_valid(self.spectra) then return
  spectra = *self.spectra
  SExcat = (*self.SExcat)[*,0]
  openw, lun, filename, /get_lun
  printf, lun, 'global color=green font="helvetica 10 normal" '+ $
          'select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source'
  printf, lun, 'fk5'
  printf, lun
  printf, lun
  for i=0, n_elements(spectra)-1 do begin
     ; use boxcircle if cluster member, otherwise use circle
     if abs(spectra[i].z-self.zcluster) lt 0.015 then begin
        pointtype='boxcircle'
     endif else begin
        pointtype='circle'
     endelse
     ; use red for spectra not found in the SExcat
     if spectra[i].galid eq -1 then begin
        printf, lun, format='(%"point(  %f15, %f15 ) # point=%s color=red")', $
                spectra[i].ra, spectra[i].dec, pointtype
     endif else begin ;use blue for spectra that were found in the SExcat
        printf, lun, format='(%"point(  %f15, %f15 ) # text={%s} point=%s color=blue")', $
                spectra[i].ra, spectra[i].dec, $
                spectra[i].clusterid+string(spectra[i].galid,format='(I05)'), $
                pointtype
        dx = -(SExcat[spectra[i].galid].alphawin_j2000*3600-spectra[i].ra*3600) $
             *cos(spectra[i].dec*!dpi/180.d)
        dy = (SExcat[spectra[i].galid].deltawin_j2000*3600-spectra[i].dec*3600)
        angle = atan(dy/dx)*180.d/!dpi
        if dx lt 0 then angle += 180.d
        length2 = dx*dx+dy*dy
        length = sqrt(length2)
        if ~finite(angle) then continue
        printf, lun, format='(%"vector( %f15, %f15, %f15\", %f15 ) # color=blue")', $
                spectra[i].ra, spectra[i].dec, length, angle
     endelse
  endfor
  close, lun
  free_lun, lun
end

Function Cluster2::Color, gals=gals, $
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
     if ~keyword_set(silent) then print, 'WARNING: using minimum radius of 3.'
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
        if n_elements(ap_factor) ne 0 then r = ap_factor*(*self.galfitcat)[gals[igal]].re $
        else r = ap_radius
        if r lt 0 then begin
           color[igal]=!values.f_nan
           continue
        endif
        r = r > min_radius
        mag1 = self->ApPhot( r, band=band1, gals=gals[igal], clean=clean, xpsf=xpsf, _extra=extra, /silent )
        mag2 = self->ApPhot( r, band=band2, gals=gals[igal], clean=clean, xpsf=xpsf, _extra=extra, /silent )
        color[igal] = mag1 - mag2
        if arg_present(mags) then begin
           mags[igal, 0] = mag1
           mags[igal, 1] = mag2
        endif
     endfor
  endif else begin
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
              re = (*self.galfitcat)[gals[igal]].re
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
              mag1 = interpol( (*self.xpsfcat)[gals[igal],iband1].cleanap, $
                               aps, ap_radius > min_radius )
              mag2 = interpol( (*self.xpsfcat)[gals[igal],iband2].cleanap, $
                               aps, ap_radius > min_radius )
              color[igal] = mag1 - mag2
              if arg_present(mags) then begin
                 mags[igal, 0] = mag1
                 mags[igal, 1] = mag2
              endif
           endelse
        endfor
     endif else begin
        for igal=0, n_elements(gals)-1 do begin
           if n_elements(ap_factor) ne 0 then radius = ap_factor*(*self.galfitcat)[gals[igal]].re $
           else radius = ap_radius
           radius = radius > min_radius
           mag1 = interpol((*self.SExcat)[gals[igal],iband1].mag_aper, diameters/2.0, radius)
           mag2 = interpol((*self.SExcat)[gals[igal],iband2].mag_aper, diameters/2.0, radius)
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

Function Cluster2::ApPhot, radius, $
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
  if arg_present(err) then errflux = fltarr(n_elements(gals),n_elements(radius))*!values.f_nan

  if keyword_set(clean) then begin
     psf = *(*self.psf)[iband]
     psf = undersampleimage( psf, self.samplefactor )
  endif else if keyword_set(xpsf) then begin
     xpsf = undersampleimage( *(*self.psf)[xband], self.samplefactor )
     xpsf /= total(xpsf)
  endif

  for igal=0, n_elements(gals)-1 do begin
     if ~keyword_set(silent) then counter, igal+1, n_elements(gals), 'Calculating aperture photometry: '
     SExcat1 = (*self.SExcat)[gals[igal],0]
     p = round([ SExcat1.xwin_image-1, $ ;;SEx coords to IDL coords
                 SExcat1.ywin_image-1, $
                 SExcat1.xwin_image-1, $
                 SExcat1.ywin_image-1 ]) + [-1,-1,1,1]*round((max(radius)*1.2) + 20)
     if total(p lt 0) ne 0 then begin
        fluxes[igal,*]=!values.d_nan
        continue
     endif
     if arg_present(err) then subimage0 = self->ImageSection( p, band=band, err=suberrimage ) $
     else subimage = self->ImageSection( p, band=band )
     sky = (*self.skycat)[gals[igal],iband].sky_value
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

Function Cluster2::CLEANCat, gals=gals, band=band, store=store
  if n_elements(gals) eq 0 then gals=lindgen(n_elements(*self.galfitcat))
  if total(gals) eq -1 then return, 0
  store=keyword_set(store)
  if ~ptr_valid(self.CLEANcat) then begin
     CLEANcat = { CLEANcat, $
                  clusterid: '', $
                  galid:     0L, $
                  band:      '', $
                  cleanap:   fltarr(5), $
                  cleanfac:  fltarr(5) }
     CLEANcat.cleanap = fltarr(5)*!values.f_nan
     CLEANcat.cleanfac = fltarr(5)*!values.f_nan
     CLEANcat = replicate(CLEANcat, n_elements(*self.galfitcat), n_elements(*self.band))
     CLEANcat.clusterid = (*self.SExcat).clusterid
     CLEANcat.galid = (*self.SExcat).galid
     self.CLEANcat = ptr_new(CLEANcat)
  endif
  CLEANcat = *self.CLEANcat
  
  if n_elements(band) eq 0 then begin
     for iband=0, n_elements(*self.band)-1 do begin
        band=(*self.band)[iband]
        CLEANcat[gals,iband] = (self->CLEANcat(gals=gals, band=band, store=store))[*,iband]
     endfor
     return, CLEANcat[gals,*]
  endif
  
  iband=(where(band eq (*self.band)))[0]
  psf = undersampleimage(*(*self.psf)[iband], self.samplefactor)
  for igal=0, n_elements(gals)-1 do begin
     counter, igal+1, n_elements(gals), 'Creating CLEANcat for band '+band+': '
     gal1 = (*self.galfitcat)[gals[igal]]
     SEx1 = (*self.SExcat)[gals[igal],0]
     sky1 = (*self.skycat)[gals[igal],iband]
     max_radius = 1.5*gal1.re > 20.
     
     p = round([ SEx1.xwin_image-1, $ ;;SEx coords to IDL coords
                 SEx1.ywin_image-1, $
                 SEx1.xwin_image-1, $
                 SEx1.ywin_image-1 ]) + [-1,-1,1,1]*round((max_radius*1.2) + 20)
     size=size(*(*self.image)[0],/dim)
     if (p[0] lt 0 or p[1] lt 0 or p[2] ge size[0] or p[3] ge size[1]) then continue
     image = self->ImageSection(p, band=band)
     if finite(sky1.sky_value) then image -= sky1.sky_value
     image = cleanimage( image, psf, 1.5 )
     apfluxes = apphot( image, $
                        SEx1.xwin_image-0.5-p[0], $ ;;SEx to josh coord
                        SEx1.ywin_image-0.5-p[1], $
                        [3., 5., 10., 15., 20.] )
     facfluxes = apphot( image, $
                         SEx1.xwin_image-0.5-p[0], $
                         SEx1.ywin_image-0.5-p[1], $
                         [0.5, 0.75, 1.0, 1.25, 1.5]*gal1.re )
     CLEANcat[gals[igal],iband].cleanap = -2.5*alog10(apfluxes)+(*self.zeropoint)[iband]
     CLEANcat[gals[igal],iband].cleanfac = -2.5*alog10(facfluxes)+(*self.zeropoint)[iband]
  endfor
  print
  if keyword_set(store) then begin
     ptr_free, self.CLEANcat
     self.CLEANcat = ptr_new(CLEANcat)
  endif
  return, CLEANcat[gals,*]
end

Function Cluster2::XPSFCat, gals=gals, band=band, cband=cband, store=store
  if n_elements(gals) eq 0 then gals=lindgen(n_elements(*self.galfitcat))
  if total(gals) eq -1 then return, 0
  store=keyword_set(store)
  if ~ptr_valid(self.XPSFcat) then begin
     XPSFcat = { XPSFcat, $
                 clusterid:  '', $
                 galid:      0L, $
                 band:       '', $
                 xpsfap:     fltarr(5), $
                 xpsffac:    fltarr(5), $
                 xpsfaperr:  fltarr(5), $
                 xpsffacerr: fltarr(5) }
     XPSFcat.xpsfap = fltarr(5)*!values.f_nan
     XPSFcat.xpsffac = fltarr(5)*!values.f_nan
     XPSFcat = replicate(XPSFcat, n_elements(*self.galfitcat), n_elements(*self.band))
     XPSFcat.clusterid = (*self.SExcat).clusterid
     XPSFcat.galid = (*self.SExcat).galid
     XPSFcat.band = (*self.SExcat).band
     self.XPSFcat = ptr_new(XPSFcat)
  endif
  XPSFcat = *self.XPSFcat
  
  if n_elements(band) eq 0 then begin
     for iband=0, n_elements(*self.band)-1 do begin
        band=(*self.band)[iband]
        XPSFcat[gals,iband] = (self->XPSFcat(gals=gals, band=band, store=store))[*,iband]
     endfor
     return, XPSFcat[gals,*]
  endif
  
  iband=(where(band eq (*self.band)))[0]
  if n_elements(cband) eq 0 then cband=(*self.band)[1-iband] ;;this works when only i & z bands present
  xband=(where(cband eq (*self.band)))[0]
  psf = undersampleimage(*(*self.psf)[iband], self.samplefactor)
  xpsf = undersampleimage(*(*self.psf)[xband], self.samplefactor)
  xpsf /= total(xpsf)
  for igal=0, n_elements(gals)-1 do begin
     counter, igal+1, n_elements(gals), 'Creating XPSFcat for band '+band+': '
     gal1 = (*self.galfitcat)[gals[igal]]
     SEx1 = (*self.SExcat)[gals[igal],0]
     sky1 = (*self.skycat)[gals[igal],iband]
     max_radius = 1.5*gal1.re > 20.
     
     p = round([ SEx1.xwin_image-1, $ ;;SEx coords to IDL coords
                 SEx1.ywin_image-1, $
                 SEx1.xwin_image-1, $
                 SEx1.ywin_image-1 ]) + [-1,-1,1,1]*round((max_radius*1.2) + 20)
     size=size(*(*self.image)[iband],/dim)
     if (p[0] lt 0 or p[1] lt 0 or p[2] ge size[0] or p[3] ge size[1]) then continue
;     image = self->ImageSection(p, band=band, err=errimage)
;     if finite(sky1.sky_value) then image -= sky1.sky_value
;     image = convol( image, xpsf, /edge_truncate )
     junk = self->ImageSection(p, band=band, err=errimage, cvlim=image)
     if finite(sky1.sky_value) then image -= sky1.sky_value
     apfluxes = apphot( image, $
                        SEx1.xwin_image-0.5-p[0], $  ;;SEx to josh coords
                        SEx1.ywin_image-0.5-p[1], $
                        [3., 5., 10., 15., 20.] )
     facfluxes = apphot( image, $
                         SEx1.xwin_image-0.5-p[0], $
                         SEx1.ywin_image-0.5-p[1], $
                         [0.5, 0.75, 1.0, 1.25, 1.5]*gal1.re )
     aperrflux = sqrt(apphot( errimage^2, $
                              SEx1.xwin_image-0.5-p[0], $  ;;SEx to josh coords
                              SEx1.ywin_image-0.5-p[1], $
                              [3., 5., 10., 15., 20.] ))
     facerrflux = sqrt(apphot( errimage^2, $
                               SEx1.xwin_image-0.5-p[0], $
                               SEx1.ywin_image-0.5-p[1], $
                               [0.5, 0.75, 1.0, 1.25, 1.5]*gal1.re ))
     XPSFcat[gals[igal],iband].xpsfap = -2.5*alog10(apfluxes)+(*self.zeropoint)[iband]
     XPSFcat[gals[igal],iband].xpsffac = -2.5*alog10(facfluxes)+(*self.zeropoint)[iband]
     XPSFcat[gals[igal],iband].xpsfaperr = $
        2.5/2.*alog10((1.+aperrflux/apfluxes)/(1.-aperrflux/apfluxes))
     XPSFcat[gals[igal],iband].xpsffacerr = $
        2.5/2.*alog10((1.+facerrflux/facfluxes)/(1.-facerrflux/facfluxes))
  endfor
  print
  if keyword_set(store) then begin
     ptr_free, self.XPSFcat
     self.XPSFcat = ptr_new(XPSFcat)
  endif
  return, XPSFcat[gals,*]
end

Function Cluster2::NearestObjs, ra, dec
  SExcat=(*self.SExcat)[*,0]
  dist2=(ra-SExcat.alphawin_j2000)^2*(cos(dec*!dpi/180.d))^2+(dec-SExcat.deltawin_j2000)^2
  s=sort(dist2)
  return, SExcat[s[0:9]].galid
end

Function Cluster2::VisMorph, filename, radec=radec
  VisMorphCat = { VisMorphCat, $
                  clusterid: '', $
                  galid:     0L, $
                  morph:     '', $
                  ttype:     99 }
  VisMorphCat = replicate(VisMorphCat, n_elements(*self.galfitcat))
  VisMorphCat.clusterid = (*self.galfitcat).clusterid
  VisMorphCat.galid = (*self.galfitcat).galid

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
        ra0 = ten(rasplit)*15.
        dec0 = ten(decsplit)
        w=self->nearestobjs(ra0, dec0)
        VisMorphcat[w[0]].morph = str(ttype[iid])
        VisMorphcat[w[0]].ttype = ttype[iid]
     endfor
  endelse
  return, VisMorphCat
end

Pro Cluster2::List, filename, gals=gals
  if n_elements(gals) eq 0 then gals = lindgen(n_elements(self.galfitcat))
  s=self->summary()
  openw, lun, filename, /get_lun
  for igal=0, n_elements(gals)-1 do begin
     iobj = gals[igal]
     radec, s[iobj].alphawin_j2000, s[iobj].deltawin_j2000, ihr, imin, xsec, ideg, imn, xsc
     rastr = string(format='(%"%02i %02i %06.3f")', ihr, imin, xsec)
     decstr = string(format='(%"%+03i %02i %06.3f")', ideg, imn, xsc)
     printf, lun, format='(%"%6s 600 %7.3f %s %s 2000.0 2000.0 0.0 0.0")', $
             self.clusterid+str(iobj), s[iobj].zmag_aper[9], rastr, decstr  ;;10 pixel radius magnitude
  endfor
  close, lun
  free_lun, lun
end

Pro Cluster2::List2, filename, gals=gals
  if n_elements(gals) eq 0 then gals = lindgen(n_elements(self.galfitcat))
  s=self->summary()
  openw, lun, filename, /get_lun
  for igal=0, n_elements(gals)-1 do begin
     iobj = gals[igal]
     printf, lun, format='(%"%6s %7.2f %7.2f %9.4f")', $
             self.clusterid+str(iobj), s[iobj].xwin_image, s[iobj].ywin_image, $
             10^(-0.4*(s[iobj].zmag_aper[9]-*(self.zeropoint[0]))) ;;10 pixel radius magnitude
  endfor
  close, lun
  free_lun, lun
end

Function Cluster2::ImportRegion, filename
  gals=[0]
  openr, lun, filename, /get_lun
  a=''
  while ~eof(lun) do begin
     readf, lun, a
     str=strsplit(a,'(,)',/extract)
     for is=0, n_elements(str)-1 do begin
        if str[0] eq 'ellipse' and strpos(str[is], 'text') ne -1 then begin
           w=strpos(str[is],'{')
           id=strmid(str[is], w+1, 5)
           gals = [gals,long(id)]
        endif
     endfor
  endwhile
  close, lun
  free_lun, lun
  return, gals[1:*]
end

Function Cluster2::wgals
  s=(*self.SExcat)[*,0]
  c=self->color(ap_radius=10)
  wgals=where( s.flux_radius gt 2.2 and s.mag_auto gt 17 and s.mag_auto lt 27 and c gt -1 and c lt 2 and s.a_image*s.kron_radius lt 200 )
  return, wgals
end

Function Cluster2::wstars
  s=(*self.SExcat)[*,0]
  wstars=where( s.flux_radius gt 1.4 and s.flux_radius lt 1.72 and s.mag_auto gt 18 and s.mag_auto lt 24 )
  return, wstars
end

Function Cluster2::sxpar, keyword, band=band
  if n_elements(band) eq 0 then band = (*self.band)[0]
  iband = (where(band eq *self.band))[0]
  hdr=headfits((*self.filehead)[iband]+'_ERR.fits.gz')
  return, sxpar( hdr, keyword )
end

Function Cluster2::hdr, band=band
  if n_elements(band) eq 0 then band = (*self.band)[0]
  iband = (where(band eq *self.band))[0]
  return, headfits((*self.filehead)[iband]+'_SCI.fits.gz')
end

Function Cluster2::Summary
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
        skycat = struct_trimtags((*self.skycat)[*,iband], except=['CLUSTERID','GALID','BAND'])
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
        CLEANcat = struct_trimtags((*self.CLEANcat)[*,iband], except=['CLUSTERID','GALID','BAND'])
        tags = tag_names(CLEANcat)
        CLEANcat = struct_renametags(CLEANcat, (*self.band)[iband]+tags)
        outcat = struct_addtags(outcat, CLEANcat)
     endfor
  endif
  if ptr_valid(self.XPSFcat) then begin
     for iband=0, nbands-1 do begin
        XPSFcat = struct_trimtags((*self.XPSFcat)[*,iband], except=['CLUSTERID','GALID','BAND'])
        tags = tag_names(XPSFcat)
        XPSFcat = struct_renametags(XPSFcat, (*self.band)[iband]+tags)
        outcat = struct_addtags(outcat, XPSFcat)
     endfor
  endif
  if ptr_valid(self.speccat) then begin
     speccat = struct_trimtags(*self.speccat, except=['CLUSTERID','GALID'])
     outcat = struct_addtags(outcat, speccat)
  endif
  return, outcat
end

Function Cluster2::Bencat
  outcat1 = { number:0L, $
              mag_auto:!values.f_nan, $
              magerr_auto:!values.f_nan, $
              background:!values.f_nan, $
              threshold:!values.f_nan, $
              mu_threshold:!values.f_nan, $
              flux_max:!values.f_nan, $
              mu_max:!values.f_nan, $
              isoarea_image:!values.f_nan, $
              isoarea_world:!values.f_nan, $
              xmin_image:0L, $
              ymin_image:0L, $
              xmax_image:0L, $
              ymax_image:0L, $
              x_image:!values.f_nan, $
              y_image:!values.f_nan, $
              x_world:!values.f_nan, $
              y_world:!values.f_nan, $
              xpeak_image:0L, $
              ypeak_image:0L, $
              xpeak_world:!values.d_nan, $
              ypeak_world:!values.d_nan, $
              alpha_j2000:!values.d_nan, $
              delta_j2000:!values.d_nan, $
              x2_image:!values.d_nan, $
              y2_image:!values.d_nan, $
              xy_image:!values.d_nan, $
              x2_world:!values.d_nan, $
              y2_world:!values.d_nan, $
              xy_world:!values.d_nan, $
              a_image:!values.f_nan, $
              b_image:!values.f_nan, $
              a_world:!values.f_nan, $
              b_world:!values.f_nan, $
              thetawin_image:!values.f_nan, $
              theta_world:!values.f_nan, $
              flags:0, $
              class_star:!values.f_nan, $
              kron_radius:!values.f_nan, $
              cxx_image:!values.f_nan, $
              cyy_image:!values.f_nan, $
              cxy_image:!values.f_nan, $
              mag_aper:!values.f_nan, $
              magerr_aper:!values.f_nan, $
              mag_iso:!values.f_nan, $
              magerr_iso:!values.f_nan, $
              z_ap:!values.f_nan, $
              i_ap:!values.f_nan, $
              iz:!values.f_nan, $
              izerr:!values.f_nan, $
              zerr_p:!values.f_nan, $
              ierr_p:!values.f_nan, $
              izerr_p:!values.f_nan, $
              dist:!values.f_nan, $
              re:!values.f_nan, $
              axrat:!values.f_nan, $
              angle:!values.f_nan, $
              sky:!values.f_nan, $
              gini:!values.f_nan, $
              asym:!values.f_nan, $
              dimgr:!values.f_nan, $
              dimgrerr:!values.f_nan }

  outcat = replicate(outcat1, n_elements(*self.galfitcat))
  sexcat = (*self.sexcat)[*,0]
  galfitcat = (*self.galfitcat)
  skycat = (*self.skycat)[*,0]
  morphcat = (*self.morphcat)
  outcat.number = sexcat.galid
  outcat.mag_auto = sexcat.mag_auto
  outcat.magerr_auto = sexcat.magerr_auto
  outcat.background = sexcat.background
  outcat.mu_threshold = sexcat.mu_threshold
  outcat.mu_max = sexcat.mu_max
  outcat.isoarea_image = sexcat.isoarea_image
  outcat.x_image = sexcat.xwin_image
  outcat.y_image = sexcat.ywin_image
  outcat.alpha_j2000 = sexcat.alphawin_j2000
  outcat.delta_j2000 = sexcat.deltawin_j2000
  outcat.a_image = sexcat.a_image
  outcat.b_image = sexcat.b_image
  outcat.thetawin_image = sexcat.theta_image
  outcat.class_star = sexcat.class_star
  outcat.kron_radius = sexcat.kron_radius
  w=where((*self.galfitcat).re ge 0.)
  color = self->color(gals=w, /xpsf, mags=mags)
  err = geterrs7(self)
  outcat[w].z_ap = mags[*,0]
  outcat[w].i_ap = mags[*,1]
  outcat[w].iz=color
  outcat.izerr=err
  outcat.re=galfitcat.re
  outcat.axrat=galfitcat.boa
  outcat.angle=galfitcat.pa
  outcat.sky=skycat.sky_value
  outcat.gini=morphcat.gini
  outcat.asym=morphcat.asym
  return, outcat
end

Pro Cluster2::CleanUp
  self->FreeImages
  ptr_free, (*self.galapagoscat).maskfit
  ptr_free, (*self.galapagoscat).simfit
  ptr_free, (*self.galapagoscat).maskfit2
  ptr_free, (*self.galapagoscat).simfit2
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
  ptr_free, self.CLEANcat
  ptr_free, self.XPSFcat
  ptr_free, self.speccat
  ptr_free, self.spectra

  ptr_free, self.image
  ptr_free, self.errim
  ptr_free, self.whtim
  ptr_free, self.skyim
  ptr_free, self.cvlim
  ptr_free, self.mskim
  ptr_free, self.hdr
  ptr_free, self.psf
end

Pro Cluster2__Define, struct
  struct = { Cluster2, $
             Inherits JEMobject, $
             clustername:   '', $
             clusterid:     '', $
             zcluster:      !values.d_nan, $
             veldisp:       !values.d_nan, $
             ra:            !values.d_nan, $
             dec:           !values.d_nan, $
             ebv:           !values.d_nan, $
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
             XPSFcat:       ptr_new(), $
             speccat:       ptr_new(), $
             spectra:       ptr_new(), $
             image:         ptr_new(), $
             errim:         ptr_new(), $
             whtim:         ptr_new(), $
             skyim:         ptr_new(), $
             cvlim:         ptr_new(), $
             mskim:         ptr_new(), $
             seghot:        ptr_new(), $
             segcold:       ptr_new(), $
             hdr:           ptr_new(), $
             psf:           ptr_new(), $
             samplefactor:  0L }
end
