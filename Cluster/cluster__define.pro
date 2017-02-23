Function Cluster::Init, $
   _extra=extra

  self->FreeImages
  self->SetProperty, _extra=extra
  return, 1
end

Pro Cluster::Sextract, $
   sexfile, $
   outdir=outdir, $
   _extra=extra

  if n_elements(outdir) eq 0 then cd, current=outdir
  outdir=directoryify(outdir)
  ;;do the sextraction
  sextract_iz, self.ifilename, self.zfilename, $
               sexfile, self.clusterid, outdir=outdir, _extra=extra
  ;;clean the results (this code is cleaner in python...)
  icat = outdir+self.clusterid+'_i.cat'
  zcat = outdir+self.clusterid+'_z.cat'
  spawn, 'cp '+icat+' '+icat+'.orig'
  spawn, 'cp '+zcat+' '+zcat+'.orig'
  spawn, 'python2.5 /home/jmeyers314/olivetos/jemidl/Cluster/cleancats.py '+ $
         icat+' '+zcat
end

Pro Cluster::SextractImport, $
   icat, zcat

  self.sexcat = ptr_new( import_sexcat( icat, zcat, self.clusterid ) )
end

Pro Cluster::SextractColdHot, $
   coldsexfile, hotsexfile, $
   outdir=outdir
  
  if n_elements(outdir) eq 0 then cd, current=outdir
  outdir=directoryify(outdir)
  ;;do the cold sextraction
  sextract_iz, self.ifilename, self.zfilename, $
               coldsexfile, self.clusterid, outdir=outdir, /segment
  ;;clean the results (this code is neater in python...)
  icat = outdir+self.clusterid+'_i.cat'
  zcat = outdir+self.clusterid+'_z.cat'
  spawn, 'python2.5 /home/jmeyers314/olivetos/jemidl/Cluster/cleancats.py '+ $
         icat+' '+zcat
  file_copy, icat, outdir+self.clusterid+'_cold_i.cat', /overwrite
  file_copy, zcat, outdir+self.clusterid+'_cold_z.cat', /overwrite
  new_icat = outdir+self.clusterid+'_cold_i.cat'
  new_zcat = outdir+self.clusterid+'_cold_z.cat'
  ;;load to tmp structure
  coldcat = import_sexcat( new_icat, new_zcat, self.clusterid )
  
  ;;now for the hot step
  sextract_iz, self.ifilename, self.zfilename, $
               hotsexfile, self.clusterid, outdir=outdir
  icat = outdir+self.clusterid+'_i.cat'
  zcat = outdir+self.clusterid+'_z.cat'
  spawn, 'python2.5 /home/jmeyers314/olivetos/jemidl/Cluster/cleancats.py '+ $
         icat+' '+zcat
  file_copy, icat, outdir+self.clusterid+'_hot_i.cat', /overwrite
  file_copy, zcat, outdir+self.clusterid+'_hot_z.cat', /overwrite
  new_icat = outdir+self.clusterid+'_hot_i.cat'
  new_zcat = outdir+self.clusterid+'_hot_z.cat'
  hotcat = import_sexcat( new_icat, new_zcat, self.clusterid )
  segment = mrdfits('/home/jmeyers314/scp2/clustersnew/'+self.clusterid+'/'+self.clusterid+'segment.fits')
  segment = segment ne 0
  for ix=-2,2 do begin
     segment = segment or shift(segment, ix, 0)
  endfor
  for iy=-2,2 do begin
     segment = segment or shift(segment, 0, iy)
  endfor
  cat = coldcat
  for ihot = 1, n_elements(hotcat)-1 do begin
     hotcat1 = hotcat[ihot]
     x=fix(hotcat1.xwin_image)
     y=fix(hotcat1.ywin_image)
     hotcat1.galid = cat[n_elements(cat)-1].galid+1
     if segment[x,y] eq 0 then cat = [cat, hotcat1]
  endfor
  self.sexcat = ptr_new( cat )
end

Pro Cluster::Galapagos
  if ptr_valid(self.galapagoscat) then ptr_free, self.galapagoscat
  self.galapagoscat = ptr_new( galapagos( *self.sexcat ) )
end

Pro Cluster::LoadImages, epoch=epoch
  self->LoadiImage, epoch=epoch
  self->LoadiErrim, epoch=epoch
  self->LoadzImage, epoch=epoch
  self->LoadzErrim, epoch=epoch
end

Function Cluster::EpochizeName, name, epoch
  p=strpos(name, 'Epoch')
  s=strmid(name, 0, p)+'Epoch'+string(epoch,format='(I02)')+strmid(name, p+7)
  return, s
end

Pro Cluster::LoadiImage, epoch=epoch
  if n_elements(epoch) eq 0 then epoch=0
  if ptr_valid(self.iimage) and self.iimepoch eq epoch then return
  ifilename = self->EpochizeName(self.ifilename, epoch)
  ptr_free, self.iimage
  self.iimage = ptr_new()
  ptr_free, self.ihead
  self.iimepoch = -1
  if ~(file_info(ifilename)).exists then return
  print, 'reading F775W image'
  iimage = mrdfits( ifilename, 1, /silent )
  ihead = headfits( ifilename, /silent )
  self.iimage = ptr_new(iimage)
  self.ihead = ptr_new(ihead)
  self.iimepoch = epoch
end

Pro Cluster::LoadzImage, epoch=epoch
  if n_elements(epoch) eq 0 then epoch=0
  if ptr_valid(self.zimage) and self.zimepoch eq epoch then return

  zfilename = self->EpochizeName(self.zfilename, epoch)
;  zfilename=self.zfilename

  ptr_free, self.zimage
  self.zimage = ptr_new()
  ptr_free, self.zhead
  self.zimepoch = -1
  if ~(file_info(zfilename)).exists then return
  print, 'reading F850LP image'
  zimage = mrdfits( zfilename, 1, /silent )
  zhead = headfits( zfilename, /silent )
  self.zimage = ptr_new(zimage)
  self.zhead = ptr_new(zhead)
  self.zimepoch = epoch
end

Pro Cluster::LoadiErrim, epoch=epoch
  if n_elements(epoch) eq 0 then epoch=0
  if ptr_valid(self.ierrim) and self.ierepoch eq epoch then return
  self->LoadiImage, epoch=epoch

  ifilename = self->EpochizeName(self.ifilename, epoch)
  ptr_free, self.ierrim
  self.ierrim = ptr_new()
  self.ierepoch = -1
  if ~(file_info(ifilename)).exists then return

  print, 'reading F775W pixel time image'
  itime = mrdfits( ifilename, 2, /silent )
  print, 'reading F775W context image'
  icont = mrdfits( ifilename, 3, /silent )

  nx = (size( *self.iimage, /dimens ))[0]
  ny = (size( *self.iimage, /dimens ))[1]

  ihead=*self.ihead

  igain = mean( double( [sxpar( ihead, 'ATODGNA' ), $
                         sxpar( ihead, 'ATODGNB' ), $
                         sxpar( ihead, 'ATODGNC' ), $
                         sxpar( ihead, 'ATODGND' )] ) )
;; read noise per exposure
  irdnoise = mean( double( [sxpar( ihead, 'READNSEA' ), $
                            sxpar( ihead, 'READNSEB' ), $
                            sxpar( ihead, 'READNSEC' ), $
                            sxpar( ihead, 'READNSED' )] ) ) 
;  isky = double( sxpar( ihead, 'SKY_MEAN' ) )  ;; counts per 500s(?)
;  isky /= 500  ;; seconds per exposure
  isky = 0.0504

  dark=0.00318 ;;counts per second

  inexp = intarr( nx, ny )
  if size( icont, /n_dim ) eq 2 then inslice=1 $
  else inslice = (size( icont, /dim ))[2]
  print, 'flattening F775W context image'
  for iislice=0, inslice-1 do begin
     for i=0, 31 do begin
        inexp += (icont[*,*,iislice] and 2L^i) eq 2L^i
     endfor
  endfor
  print, 'computing F775W error map'
  ierrim = sqrt( *self.iimage*itime*igain $
                 + irdnoise*irdnoise*inexp $
                 + (isky + dark)*itime*igain ) / (igain*itime)
  ptr_free, self.ierrim
  self.ierrim = ptr_new(ierrim)
  self.ierepoch = epoch
end

Pro Cluster::LoadzErrim, epoch=epoch
  if n_elements(epoch) eq 0 then epoch=0
  if ptr_valid(self.zerrim) and self.zerepoch eq epoch then return
  self->LoadzImage, epoch=epoch

  zfilename = self->EpochizeName(self.zfilename, epoch)
;  zfilename = self.zfilename

  ptr_free, self.zerrim
  self.zerrim = ptr_new()
  self.zerepoch = -1
  if ~(file_info(zfilename)).exists then return

  print, 'reading F850LP pixel time image'
  ztime = mrdfits( zfilename, 2, /silent )
  print, 'reading F850LP context image'
  zcont = mrdfits( zfilename, 3, /silent )

  nx = (size( *self.zimage, /dimens ))[0]
  ny = (size( *self.zimage, /dimens ))[1]

  zhead = *self.zhead

  zgain = mean( double( [sxpar( zhead, 'ATODGNA' ), $
                         sxpar( zhead, 'ATODGNB' ), $
                         sxpar( zhead, 'ATODGNC' ), $
                         sxpar( zhead, 'ATODGND' )] ) )
;; read noise per exposure
  zrdnoise = mean( double( [sxpar( zhead, 'READNSEA' ), $
                            sxpar( zhead, 'READNSEB' ), $
                            sxpar( zhead, 'READNSEC' ), $
                            sxpar( zhead, 'READNSED' )] ) )
;  zsky = double( sxpar( zhead, 'SKY_MEAN' ) )  ;; counts per 500s(?)
;  zsky /= 500  ;; seconds per exposure
  zsky = 0.0268

  dark=0.00318 ;; counts per second

  znexp = intarr( nx, ny )
  if size( zcont, /n_dim ) eq 2 then znslice=1 $
  else znslice = (size( zcont, /dim ))[2]
  print, 'flattening F850LP context image'
  for izslice=0, znslice-1 do begin
     for i=0, 31 do begin
        znexp += (zcont[*,*,izslice] and 2L^i) eq 2L^i
     endfor
  endfor
  print, 'computing F850LP error map'
  zerrim = sqrt( *self.zimage*ztime*zgain $
                 + zrdnoise*zrdnoise*znexp $
                 + (zsky + dark)*ztime*zgain ) / (zgain*ztime)
  ptr_free, self.zerrim
  self.zerrim = ptr_new(zerrim)
  self.zerepoch = epoch
end

Function Cluster::ImageSection, $
   p, $      ;; position to extract
   band=band, $
   err=err, $
   psfmatch=psfmatch, $
   CLEAN=CLEAN

  if n_elements(band) eq 0 then band='z'
  case band of 
     'z' : begin
        if ~ptr_valid(self.zimage) then begin
           out=dblarr(p[2]-p[0]+1,p[3]-p[1]+1)*!values.d_nan
        endif else begin
           if arg_present(err) then begin
              if ~ptr_valid(self.zerrim) then $
                 err=dblarr(p[2]-p[0]+1,p[3]-p[1]+1)*!values.d_nan $
              else $
                 err = (*self.zerrim)[p[0]:p[2],p[1]:p[3]]
           endif
           out=(*self.zimage)[p[0]:p[2],p[1]:p[3]]
        endelse
     end
     'i' : begin
        if ~ptr_valid(self.iimage) then begin
           out=dblarr(p[2]-p[0]+1,p[3]-p[1]+1)*!values.d_nan
        endif else begin
           if arg_present(err) then begin
              if ~ptr_valid(self.ierrim) then $
                 err=dblarr(p[2]-p[0]+1,p[3]-p[1]+1)*!values.d_nan $
              else $
                 err = (*self.ierrim)[p[0]:p[2],p[1]:p[3]]
           endif
           out=(*self.iimage)[p[0]:p[2],p[1]:p[3]]
        endelse
     end
  endcase

  if keyword_set(psfmatch) then begin
     case band of 
        'i' : opppsf = undersampleimage(*self.zpsf,self.psffactor)
        'z' : opppsf = undersampleimage(*self.ipsf,self.psffactor)
     endcase
     opppsf /= max(opppsf)
     out = convol(out,opppsf,total(opppsf),/edge_zero)
  endif
  if keyword_set(CLEAN) then begin
     case band of
        'i' : selfpsf = undersampleimage(*self.ipsf,self.psffactor)
        'z' : selfpsf = undersampleimage(*self.zpsf,self.psffactor)
     endcase
     selfpsf /= max(selfpsf)
     out=clean(out,selfpsf,1.5)
  endif
  return, out
end

Function Cluster::GalImage, gal, xsize=xsize, ysize=ysize, _extra=extra
  if n_elements(xsize) eq 0 or n_elements(ysize) eq 0 then begin
     galapcat=(*self.galapagoscat)[gal]
     p=[galapcat.xmin,galapcat.ymin,galapcat.xmax,galapcat.ymax]
  endif else begin
     s=(*self.sexcat)[gal]
     x=s.xwin_image
     y=s.ywin_image
     p=[x-xsize/2, y-ysize/2, x+xsize/2, y+ysize/2]
  endelse
  return, self->ImageSection(p,_extra=extra)
end

Pro Cluster::FreeImages
  ptr_free, self.iimage
  ptr_free, self.zimage
  ptr_free, self.ierrim
  ptr_free, self.zerrim
  self.iimepoch=-1
  self.zimepoch=-1
  self.ierepoch=-1
  self.zerepoch=-1
end

Pro Cluster::Isophot_Sky, gals=gals, band=band
  platescale = sxpar(*self.zhead, 'D001SCAL')*0.05
  galapagoscat = *self.galapagoscat
  if n_elements(gals) eq 0 then gals=indgen(n_elements(galapagoscat))
  print, 'Isophot_Sky on band '+band
  for igal=0, n_elements(gals)-1 do begin
     print, igal+1, ' of ', n_elements(gals)
     case band of 
        'i' : begin
           gal_isophot_sky, self.iimage, self.ierrim, $
                            galapagoscat[gals[igal]], galapagoscat, $
                            platescale=platescale, $
                            sky_value=sky_value, sky_sigma=sky_sigma, sky_median=sky_median
           galapagoscat[gals[igal]].isky=sky_value
           galapagoscat[gals[igal]].isky_sigma=sky_sigma
        end
        'z' : begin
           gal_isophot_sky, self.zimage, self.zerrim, $
                            galapagoscat[gals[igal]], galapagoscat, $
                            platescale=platescale, $
                            sky_value=sky_value, sky_sigma=sky_sigma, sky_median=sky_median
           galapagoscat[gals[igal]].zsky=sky_value
           galapagoscat[gals[igal]].zsky_sigma=sky_sigma
        end
     endcase
  endfor
  ptr_free, self.galapagoscat
  self.galapagoscat = ptr_new(galapagoscat)
end

Function Cluster::GetStars, band=band, starsize=starsize
  if n_elements(band) eq 0 then band='z'
  sexcat = *self.sexcat
  platefrac = sxpar(*self.zhead, 'D001SCAL')
  platescale = platefrac*0.05
  if n_elements(starsize) eq 0 then begin
     starsize=fix(41/platefrac)
     if fix(starsize/2.) eq starsize/2. then starsize+=1 ;;ensure odd
  endif
  innersize = fix(starsize*0.23)
  if fix(innersize/2.) eq innersize/2. then innersize+=1 ;;ensure odd

  w = where( sexcat.zflux_radius gt 1.4/platefrac and $
             sexcat.zflux_radius lt 1.72/platefrac and $
             sexcat.zmag_auto gt 18 and $
             sexcat.zmag_auto lt 24 )
  halfsize=(starsize-1)/2
  onestar = {image:dblarr(starsize,starsize), $
             model:dblarr(innersize,innersize), $
             errim:dblarr(starsize,starsize), $
             modelparams: dblarr(7), $
             x_sex:0.d, $
             y_sex:0.d }
  for iw=0, n_elements(w)-1 do begin
     star = onestar
     p = round([sexcat[w[iw]].xwin_image-1, $ ;;need to take sex coords -> IDL coords
                sexcat[w[iw]].ywin_image-1, $
                sexcat[w[iw]].xwin_image-1, $
                sexcat[w[iw]].ywin_image-1])+[-1,-1,1,1]*halfsize
     star.image = self->ImageSection(p, band=band, err=err)
     if n_elements(err) eq 0 || total(~finite(err)) ne 0 then continue
     case band of
        'i' : star.image -= sexcat[w[iw]].ibackground
        'z' : star.image -= sexcat[w[iw]].zbackground
     endcase
     star.errim = err
     star.x_sex = sexcat[w[iw]].xwin_image-fix(sexcat[w[iw]].xwin_image)
     star.y_sex = sexcat[w[iw]].ywin_image-fix(sexcat[w[iw]].ywin_image)

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

Pro Cluster::Galfit, $
   gals=gals, $
   outdir=outdir, $
   _extra=extra
;; extra includes:
;; noerr=noerr
;; take2=take2 for 2nd chance GALFITing.

  if n_elements(outdir) eq 0 then cd, current=outdir
  if n_elements(gals) eq 0 then gals = lindgen( n_elements( *self.galapagoscat ) )

  if ~ptr_valid(self.galfitcat) then begin
     galfitcat1 = { galfitcat, $
                    clusterid:  '', $
                    galid:      0L, $
                    xpos:       -1.d, $
                    ypos:       -1.d, $
                    mag:        -1.d, $
                    R_e:        -1.d, $
                    boa:        -1.d, $
                    pa:         -1.d, $
                    n:          -1.d, $
                    errxpos:    -1.d, $
                    errypos:    -1.d, $
                    errmag:     -1.d, $
                    errR_e:     -1.d, $
                    errboa:     -1.d, $
                    errpa:      -1.d, $
                    errn:       -1.d, $
                    chi2:       -1.d, $
                    dof:        -1, $
                    chi2perdof: -1.d }
     galfitcat = replicate( galfitcat1, n_elements(*self.galapagoscat) )
     galfitcat.clusterid = (*self.galapagoscat).clusterid
     galfitcat.galid = (*self.galapagoscat).galid
     self.galfitcat = ptr_new(galfitcat)
  endif

  galfitcat = (*self.galfitcat)
  zpsf = undersampleimage(*self.zpsf, self.psffactor)
  zpsf /= max(zpsf)
  psf_ptr = ptr_new(zpsf)
  for igal=0L, n_elements(gals)-1 do begin
     galap1 = (*self.galapagoscat)[gals[igal]]
     mkhdr, hdr, 4, [1,1]
     sxaddpar, hdr, 'EXPTIME', sxpar(*self.zhead, 'EXPTIME')
     sxaddpar, hdr, 'GAIN', 1.0
     sxaddpar, hdr, 'RDNOISE', 5.0
     sxaddpar, hdr, 'NCOMBINE', sxpar(*self.zhead, 'NDRIZIM')/2
     galfit1 = MakeGalfit( self.zimage, self.zerrim, psf_ptr, $
                           galap1, *self.galapagoscat, 'z', $
                           outdir=outdir, header=hdr, _extra=extra )
     galfit1->Execute
     results = galfit1->Results()
     if size(results, /tname) eq 'STRUCT' then begin
        results.xpos += galap1.xmin
        results.ypos += galap1.ymin
        dest = galfitcat[gals[igal]]
        struct_assign, results, dest, /nozero
        galfitcat[gals[igal]] = dest
     endif
     obj_destroy, galfit1
  endfor
  ptr_free, psf_ptr
  ptr_free, self.galfitcat
  self.galfitcat = ptr_new(galfitcat)
end

Function Cluster::GalfitObj, $
   gal, $
   band=band, $
   outdir=outdir, $
   _extra=extra
;; extra includes:
;; noerr=noerr
;; take2=take2 for 2nd chance GALFITing.

  if n_elements(outdir) eq 0 then cd, current=outdir
  case band of 
     'i': begin 
        psf = undersampleimage(*self.ipsf, self.psffactor)
        header = *self.ihead
        imageptr = self.iimage
        errorptr = self.ierrim
     end
     'z': begin 
        psf = undersampleimage(*self.zpsf, self.psffactor)
        header = *self.zhead
        imageptr = self.zimage
        errorptr = self.zerrim
     end
  endcase
  psf /= max(psf)
  psf_ptr = ptr_new(psf)
  mkhdr, hdr, 4, [1,1]
  galap1 = (*self.galapagoscat)[gal]
  sxaddpar, hdr, 'EXPTIME', sxpar(header, 'EXPTIME')
  sxaddpar, hdr, 'GAIN', 1.0
  sxaddpar, hdr, 'RDNOISE', 5.0
  sxaddpar, hdr, 'NCOMBINE', sxpar(header, 'NDRIZIM')/2
  galfit1 = MakeGalfit( imageptr, errorptr, psf_ptr, $
                        galap1, *self.galapagoscat, band, $
                        outdir=outdir, noerr=noerr, header=hdr )
  return, galfit1
end


;this procedure will populate the skycat catalog
Pro Cluster::EpochSkyCat, gals=gals
  if n_elements(gals) eq 0 then gals = lindgen(n_elements(*self.sexcat))
  self->CalcMaxEpoch
  if ~ptr_valid(self.skycat) then begin
     skycat1 = { clusterid:(*self.sexcat)[0].clusterid, $
                 galid: 0, $
                 isky: dblarr(self.maxepoch), $
                 zsky: dblarr(self.maxepoch) }
     skycat = replicate(skycat1, n_elements(*self.sexcat))
     skycat.galid = lindgen(n_elements(*self.sexcat))
  endif else skycat = *self.skycat
  for igal=0,n_elements(gals)-1 do begin
     skycat[gals[igal]].isky *= !values.d_nan
     skycat[gals[igal]].zsky *= !values.d_nan
  endfor
  for iepoch=1, self.maxepoch do begin
     self->LoadImages, epoch=iepoch ;;if fail then (i/z)imepoch=-1
     if self.iimepoch ne -1 then begin
        for igal=0,n_elements(gals)-1 do begin
           gal_isophot_sky, self.iimage, self.ierrim, $
                            (*self.galapagoscat)[gals[igal]], *self.galapagoscat, $
                            sky_value=sky_value
           skycat[gals[igal]].isky[iepoch-1]=sky_value
        endfor
     endif
     if self.zimepoch ne -1 then begin
        for igal=0,n_elements(gals)-1 do begin
           gal_isophot_sky, self.zimage, self.zerrim, $
                            (*self.galapagoscat)[gals[igal]], *self.galapagoscat, $
                            sky_value=sky_value
           skycat[gals[igal]].zsky[iepoch-1]=sky_value
        endfor
     endif
  endfor
  ptr_free, self.skycat
  self.skycat = ptr_new(skycat)
end

Pro Cluster::Morphology, $
   gals=gals, $
   outdir=outdir
  if ~ptr_valid(self.morphcat) then begin
     morphcat = { morphcat, $
                  clusterid: '', $
                  galid:     0L, $
                  asym:      -1.d, $
                  conc:      -1.d, $
                  gini:      -1.d, $
                  radius:    -1.d, $
                  thresh:    -1.d, $
                  resid:     -1.d, $
                  bright:    -1.d, $
                  npix:      -1L }
     morphcat = replicate(morphcat, n_elements(*self.galfitcat))
     morphcat.clusterid = (*self.galfitcat).clusterid
     morphcat.galid = (*self.galfitcat).galid
     self.morphcat = ptr_new(morphcat)
  endif
  morphcat = *self.morphcat
  s = self->summary()

  if n_elements(outdir) eq 0 then cd, current=outdir
  outdir=directoryify(outdir)
  if n_elements(gals) eq 0 then gals = lindgen( n_elements( *self.galfitcat ) )
  for igal=0, n_elements(gals)-1 do begin
     s1=s[gals[igal]]
     prefix=self.clusterid+string(s1.galid,format='(I05)')+'z'
     print, prefix
     file=outdir+prefix+'.fits'
     if ~(file_info(file)).exists then continue
     struct=mrdfits(outdir+prefix+'.fits',1,/silent)
     image = struct.intsb/sxpar(*self.zhead,'EXPTIME') ;; recall intsb is already sky-subtracted
     thresh=s1.zsky_sigma

     mask = image*(1-struct.mskim) ge thresh*1.5
     npix = total((image)*(1-struct.mskim) ge thresh*1.5)
     radius = sqrt( npix / struct.boa / !dpi )

     nx=(size(image, /dim))[0]
     ny=(size(image, /dim))[1]
     dist_ellipse, ellipse, [nx, ny], $
                   struct.xpos-1, struct.ypos-1, $
                   1./struct.boa, struct.pa
     outer_ellipse = ellipse lt radius
     npix = total(outer_ellipse*mask)
     if npix eq 0 then begin
        morphcat[gals[igal]].radius=radius
        morphcat[gals[igal]].thresh=thresh
        continue
     endif
     outer_total = double(total(image*outer_ellipse*mask))
     ;;concentration
     inner_ellipse = ellipse lt 0.3*radius
     if total(inner_ellipse) eq 0 then conc = -1.d $
     else conc = total(image*inner_ellipse*mask)/outer_total
     ;;asymmetry
     rotimage=rot(image, 180., 1., $
                  struct.xpos-1, struct.ypos-1, cubic=0.8, /pivot)
     resid = 0.5d*total(abs(rotimage-image)*outer_ellipse*mask)
     correction = s1.zsky_sigma*npix/sqrt(!dpi)
     asym = (resid-correction)/outer_total
     ;;gini
     list = (image)[where(outer_ellipse and mask)]
     list = list[sort(list)]
     if npix gt 2 then begin
        gini = total(((dindgen(npix)+1)*2.-npix-1)*list)
        gini /= mean(list)*npix*(npix-1)
     endif else gini=-1.d
     morphcat[gals[igal]].asym=asym
     morphcat[gals[igal]].conc=conc
     morphcat[gals[igal]].gini=gini
     morphcat[gals[igal]].radius=radius
     morphcat[gals[igal]].thresh=thresh
     morphcat[gals[igal]].resid=resid
     morphcat[gals[igal]].bright=outer_total
     morphcat[gals[igal]].npix=npix
  endfor
  ptr_free, self.morphcat
  self.morphcat = ptr_new(morphcat)
end

Pro Cluster::AddSpectroscopy, $
   dir  ;; directory to add
  spectra = addspectra( dir, self.sexcat )
  if ptr_valid(self.spectra) then begin
     oldspectra=*self.spectra
     spectra = [oldspectra, spectra]
     ptr_free, self.spectra
     self.spectra = ptr_new(spectra)
  endif else self.spectra = ptr_new(spectra)
  self->UpdateSpecCat
end

Pro Cluster::UpdateSpecCat
  speccat = { SpecCat, $
              galid:0L, $
              z:-2.d, $
              zerr:-1.d, $
              zqual:'N', $
              comment:'' }
  speccat = replicate( speccat, n_elements(*self.sexcat) )
  speccat.galid=(*self.sexcat).galid
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

Pro Cluster::CreateRegion, filename, w=w
  if n_elements(w) eq 0 then w=lindgen( n_elements(*self.sexcat) )
  openw, lun, filename, /get_lun
  printf, lun, 'global color=green font="helvetica 10 normal" '+ $
          'select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source'
  printf, lun, 'fk5'
  printf, lun
  printf, lun
  for i=0, n_elements(w)-1 do begin
     scat = (*self.sexcat)[w[i]]
     gcat = (*self.galapagoscat)[w[i]]
     printf, lun, 'ellipse(' $
             +string(scat.alphawin_j2000,format='(D11.5)')+',' $
             +string(scat.deltawin_j2000,format='(D11.5)')+',' $
             +string(scat.a_image*scat.kron_radius*0.05,format='(D8.2)')+'",' $
             +string(scat.b_image*scat.kron_radius*0.05,format='(D8.2)') $
             +'", '+string(scat.theta_image,format='(D8.2)') $
             +') # tag={aperture} text={' $
             +string(scat.galid,format='(I05)')+'} color=yellow'
     printf, lun
  endfor
  close, lun
  free_lun, lun
end

Pro Cluster::SpecRegionCheck, filename
  if ~ptr_valid(self.spectra) then return
  spectra = *self.spectra
  sexcat = *self.sexcat
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
     ; use red for spectra not found in the sexcat
     if spectra[i].galid eq -1 then begin
        printf, lun, format='(%"point(  %f15, %f15 ) # point=%s color=red")', spectra[i].ra, spectra[i].dec, pointtype
     endif else begin ;use blue for spectra that were found in the sexcat
        printf, lun, format='(%"point(  %f15, %f15 ) # text={%s} point=%s color=blue")', $
                spectra[i].ra, spectra[i].dec, $
                spectra[i].clusterid+string(spectra[i].galid,format='(I05)'), $
                pointtype
        dx = -(sexcat[spectra[i].galid].alphawin_j2000*3600-spectra[i].ra*3600)*cos(spectra[i].dec*!dpi/180.d)
        dy = (sexcat[spectra[i].galid].deltawin_j2000*3600-spectra[i].dec*3600)
        angle = atan(dy/dx)*180.d/!dpi
        if dx lt 0 then angle += 180.d
        length2 = dx*dx+dy*dy
        length = sqrt(length2)
        if ~finite(angle) then continue
        printf, lun, format='(%"vector( %f15, %f15, %f15\", %f15 ) # color=blue")', spectra[i].ra, spectra[i].dec, length, angle
     endelse
  endfor
  close, lun
  free_lun, lun
end

Function Cluster::Color, ap_radius=ap_radius, $
                         ap_factor=ap_factor, $
                         psfmatch=psfmatch, $
                         joshmag=joshmag, $
                         clean=clean, $
                         clncat=clncat, $
                         gals=gals, $
                         _extra=extra

  if n_elements(gals) eq 0 then gals=indgen(n_elements(*self.sexcat))
  color=dblarr( n_elements(gals) )
  if n_elements(ap_factor) eq 0 and n_elements(ap_radius) eq 0 then ap_factor=1.
  diameters=[1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0, $
             11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0, $
             22.0,24.0,26.0,28.0,30.0,32.0,40.0,60.0,80.0,120.0,160.0,220.0]

  psfmatch = keyword_set(psfmatch)
  clean = keyword_set(clean)
  if keyword_set(psfmatch) or keyword_set(clean) then joshmag=1
  if keyword_set(joshmag) then begin
     iprint=0
     for igal=0, n_elements(gals)-1 do begin
        if n_elements(ap_radius) ne 0 then r=ap_radius $
        else r = ap_factor*(*self.galfitcat)[gals[igal]].r_e
        if r lt 0 then begin
           color[igal]=!values.f_nan
           continue
        endif
        if igal*100L/n_elements(gals) gt iprint then begin
           print, string(format='(%"%i%%")',iprint)
           iprint += 1
        endif
        imag = self->ApPhot(r, band='i',gals=gals[igal], psfmatch=psfmatch, clean=clean, _extra=extra)
        zmag = self->ApPhot(r, band='z',gals=gals[igal], psfmatch=psfmatch, clean=clean, _extra=extra)
        color[igal]= imag-zmag
     endfor
     return, color
  endif

  if keyword_set(clncat) then begin
     for igal=0, n_elements(gals)-1 do begin
        if n_elements(ap_factor) ne 0 then begin
           apfacs = [0.5,0.75,1.0,1.25,1.5]
           color[igal] = interpol( (*self.cleancat)[gals[igal]].cleanfac, $
                                   apfacs, ap_factor)
        endif else begin
           aps = [3., 5., 10., 15., 20.]
           color[igal] = interpol( (*self.cleancat)[gals[igal]].cleanap, $
                                   aps, ap_radius)
        endelse
     endfor
     return, color
  endif

  if n_elements(ap_factor) ne 0 then begin
     for igal=0, n_elements(gals)-1 do begin
        d = 2.*ap_factor*(*self.galfitcat)[gals[igal]].r_e
        if d lt 0 then begin
           color[igal]=!values.f_nan
           continue
        endif
        imag = interpol( (*self.sexcat)[gals[igal]].imag_aper, $
                         diameters, d )
        zmag = interpol( (*self.sexcat)[gals[igal]].zmag_aper, $
                         diameters, d )
        color[igal] = imag - zmag
     endfor
  endif else begin
     for igal=0, n_elements(gals)-1 do begin
        imag = interpol( (*self.sexcat)[gals[igal]].imag_aper, diameters, 2.*ap_radius )
        zmag = interpol( (*self.sexcat)[gals[igal]].zmag_aper, diameters, 2.*ap_radius )
        color[igal] = imag - zmag
     endfor
  endelse
  return, color
end

Pro Cluster::CalcMaxEpoch, force=force
  if self.maxepoch ne 0 and ~keyword_set(force) then return
  p=strpos(self.ifilename,'/',/reverse_search)
  dir = strmid(self.ifilename, 0, p+1)
  ifiles = file_search(dir+"*v45*Epoch*775*fits")
  zfiles = file_search(dir+"*v45*Epoch*850*fits")
  zfiles = zfiles[sort(zfiles)]
  ifiles = ifiles[sort(ifiles)]
  nz = n_elements(zfiles)
  ni = n_elements(ifiles)
  p=strpos(zfiles[nz-1],'Epoch')
  zmaxepoch = fix(strmid(zfiles[nz-1], p+5, 2))
  p=strpos(ifiles[ni-1],'Epoch')
  imaxepoch = fix(strmid(ifiles[ni-1], p+5, 2))
  maxepoch = max([imaxepoch,zmaxepoch])
  self.maxepoch = maxepoch
end

Function Cluster::LightCurves, band=band, gals=gals, radius=radius, _extra=extra
  if n_elements(gals) eq 0 then gals=lindgen(n_elements(*self.sexcat))
  
  self->CalcMaxEpoch
  lightcurves = dblarr(n_elements(gals), n_elements(radius), self.maxepoch+1)
  for iepoch=0, self.maxepoch do begin
     case band of 
        'i' : begin
           self->LoadiImage, epoch=iepoch
           self->LoadzErrIm, epoch=iepoch
        end
        'z' : begin
           self->LoadzImage, epoch=iepoch
           self->LoadzErrIm, epoch=iepoch
        end
     endcase
     self->isophot_sky, gals
     lightcurves[*,*,iepoch] = self->ApPhot(radius, band=band, gals=gals, /skysub, _extra=extra)
  endfor
  return, lightcurves
end

Pro Cluster::CalcLCcat, gals=gals, ap_radius=ap_radius, ap_factor=ap_factor
  if (n_elements(ap_radius) eq 0 and n_elements(ap_factor) eq 0) then ap_factor=1.
  LCcat1 = { clusterid:     (*self.sexcat)[0].clusterid, $
             galid:         0L, $
             zmag:          dblarr(self.maxepoch), $
             imag:          dblarr(self.maxepoch), $
             zmag_PSFmatch: dblarr(self.maxepoch), $
             imag_PSFmatch: dblarr(self.maxepoch), $
             zmag_CLEAN:    dblarr(self.maxepoch), $
             imag_CLEAN:    dblarr(self.maxepoch) }
  LCcat = replicate(LCcat1,n_elements(gals))
  LCcat.galid=gals
  for iepoch=1,self.maxepoch do begin
     self->LoadImages, epoch=iepoch
     for igal=0, n_elements(gals)-1 do begin
        if n_elements(ap_factor) ne 0 then radius=(*self.galfitcat)[gals[igal]].r_e > 3. $
        else radius=ap_radius
        LCcat[igal].zmag[iepoch-1]=self->ApPhot(radius, band='z', gals=[gals[igal]], /skysub)
        LCcat[igal].zmag_PSFmatch[iepoch-1]=self->ApPhot(radius, band='z', gals=[gals[igal]], /skysub, /PSFmatch)
        LCcat[igal].zmag_CLEAN[iepoch-1]=self->ApPhot(radius, band='z', gals=[gals[igal]], /skysub, /clean)
        LCcat[igal].imag[iepoch-1]=self->ApPhot(radius, band='i', gals=[gals[igal]], /skysub)
        LCcat[igal].imag_PSFmatch[iepoch-1]=self->ApPhot(radius, band='i', gals=[gals[igal]], /skysub, /PSFmatch)
        LCcat[igal].imag_CLEAN[iepoch-1]=self->ApPhot(radius, band='i', gals=[gals[igal]], /skysub, /clean)
     endfor
  endfor
  ptr_free, self.LCcat
  self.LCcat = ptr_new(LCcat)
end

Function Cluster::GetGalImages, gal, radius=radius, _extra=extra
  if n_elements(radius) eq 0 then radius=10
  sexcat1 = (*self.sexcat)[gal]
  p = round([ sexcat1.xwin_image-1, $ ;;sex coords to IDL coords
              sexcat1.ywin_image-1, $
              sexcat1.xwin_image-1, $
              sexcat1.ywin_image-1 ]) + [-1,-1,1,1]*round((max(radius)*1.2) + 20)
  self->CalcMaxEpoch
  for i=0, self.maxepoch do begin
     self->LoadzImage, epoch=i
     im = self->ImageSection(p,band='z',_extra=extra)
     if n_elements(ims) eq 0 then ims=im else ims = [[[ims]],[[im]]]
  endfor
  return, ims
end

Function Cluster::ApPhot, radius, $
                          band=band, $
                          gals=gals, $
                          skysub=skysub, $
                          sexsub=sexsub, $
                          _extra=extra

  if keyword_set(sexsub) then skysub=1
  if n_elements(gals) eq 0 then gals=lindgen(n_elements(*self.sexcat))
  fluxes = dblarr(n_elements(gals),n_elements(radius))
  if (band eq 'i' and self.iimepoch eq -1) or (band eq 'z' and self.zimepoch eq -1) then begin
     return, fluxes * !values.d_nan
  endif
  for igal=0, n_elements(gals)-1 do begin
     sexcat1 = (*self.sexcat)[gals[igal]]
     p = round([ sexcat1.xwin_image-1, $ ;;sex coords to IDL coords
                 sexcat1.ywin_image-1, $
                 sexcat1.xwin_image-1, $
                 sexcat1.ywin_image-1 ]) + [-1,-1,1,1]*round((max(radius)*1.2) + 20)
     if total(p lt 0) ne 0 then begin
        fluxes[igal,*]=!values.d_nan
        continue
     endif
     subimage = self->ImageSection(p, band=band, _extra=extra)
     fluxes[igal,*]=apphot(subimage, $
                           (sexcat1.xwin_image-0.5)-p[0], $
                           (sexcat1.ywin_image-0.5)-p[1], $
                           radius)
     if keyword_set(skysub) then begin
        case band of
           'i': begin
              if self.iimepoch eq 0 then begin
                 if keyword_set(sexsub) then sky = (*self.sexcat)[gals[igal]].ibackground $
                 else sky = ((*self.galapagoscat)[gals[igal]]).isky
              endif else $
                 sky = ((*self.skycat)[gals[igal]]).isky[self.iimepoch-1]
           end
           'z': begin
              if self.zimepoch eq 0 then begin
                 if keyword_set(sexsub) then sky = (*self.sexcat)[gals[igal]].zbackground $
                 else sky = ((*self.galapagoscat)[gals[igal]]).zsky
              endif else $
                 sky = ((*self.skycat)[gals[igal]]).zsky[self.zimepoch-1]
           end
        endcase
        if sky eq -1 then sky=0.d
        sky *= radius*radius*!dpi
        fluxes[igal,*] -= sky
     endif
  endfor
  case band of 
     'i': begin
        mags = -2.5*alog10(fluxes)+25.678
     end
     'z': begin
        mags = -2.5*alog10(fluxes)+24.867 
     end
  endcase
  return, mags
end

Function Cluster::ApEllipPhot, radius, band=band, gals=gals, $
                               _extra=extra
  if n_elements(gals) eq 0 then gals=lindgen(n_elements(*self.sexcat))
  fluxes = dblarr(n_elements(gals),n_elements(radius))
  for igal=0,n_elements(gals)-1 do begin
     galfitcat1 = (*self.galfitcat)[gals[igal]]
     p = fix([ galfitcat1.zxpos-0.5, $
               galfitcat1.zypos-0.5, $
               galfitcat1.zxpos-0.5, $
               galfitcat1.zypos-0.5 ]) + [-1,-1,1,1]*fix((max(radius)*1.2) + 20)
     if total(p lt 0) ne 0 then begin
        fluxes[igal,*]=dblarr(n_elements(radius))*(-1.)
        continue
     endif
     subimage = self->ImageSection(p, band=band, _extra=extra)
     fluxes[igal,*]=apellipphot(subimage, $
                                (galfitcat1.zxpos-0.5)-p[0], $
                                (galfitcat1.zypos-0.5)-p[1], $
                                galfitcat1.zboa, $
                                galfitcat1.zpa+90, $
                                radius)
  endfor
  return, fluxes
end

Function Cluster::PhotDiag, gal
  if n_elements(band) eq 0 then band='z'
  self->CalcMaxEpoch
  for i=1, self.maxepoch do begin
     self->LoadzImage, epoch=i
     self->LoadzErrim, epoch=i
     gal_isophot_sky, self.zimage, self.zerrim, (*self.galapagoscat)[gal], *self.galapagoscat, $
                      sky_value=sky_value, $
                      sky_sigma=sky_sigma, $
                      sky_median=sky_median, $
                      skyim=skyim, $
                      skyreg=skyreg
     skystruct={ skyim:skyim, $
                 skyreg:skyreg, $
                 sky_value:sky_value, $
                 sky_sigma:sky_sigma, $
                 sky_median:sky_median }
     if n_elements(skystructs) eq 0 then begin
        skystructs=skystruct
     endif else begin
        skystructs=[skystructs, skystruct]
     endelse
  endfor
  return, skystructs
end

Function Cluster::Summary
  struct = { zcluster:self.zcluster, $
             veldisp:self.veldisp, $
             velocity: 0.d }
  struct = replicate(struct, n_elements(*self.sexcat))
  outcat = struct_addtags(*self.sexcat, struct)
  if ptr_valid(self.galapagoscat) then outcat = struct_addtags(outcat, struct_trimtags(*self.galapagoscat, except=['CLUSTERID','GALID']))
  if ptr_valid(self.galfitcat) then outcat = struct_addtags(outcat, struct_trimtags(*self.galfitcat, except=['CLUSTERID','GALID']))
  if ptr_valid(self.morphcat) then outcat = struct_addtags(outcat, struct_trimtags(*self.morphcat, except=['CLUSTERID','GALID']))
  self->UpdateSpecCat
  if ptr_valid(self.speccat) then outcat = struct_addtags(outcat, struct_trimtags(*self.speccat, except=['CLUSTERID','GALID']))
  if ptr_valid(self.CLEANcat) then outcat = struct_addtags(outcat, struct_trimtags(*self.CLEANcat, except=['CLUSTERID','GALID']))
  return, outcat
end

Function Cluster::NearestObjs, ra, dec
  sexcat=*self.sexcat
  dist2=(ra-sexcat.alphawin_j2000)^2*(cos(dec*!dpi/180.d))^2+(dec-sexcat.deltawin_j2000)^2
  s=sort(dist2)
  return, sexcat[s[0:9]].galid
end

Pro Cluster::CleanUp
  self->FreeImages
  ptr_free, self.sexcat
  ptr_free, (*self.galapagoscat).maskfit
  ptr_free, (*self.galapagoscat).simfit
  ptr_free, self.galapagoscat
  ptr_free, self.galfitcat
  if ptr_valid(self.spectra) then obj_destroy, (*self.spectra).spectrum
  ptr_free, self.morphcat
  ptr_free, self.CLEANcat
  ptr_free, self.spectra
  ptr_free, self.speccat
  ptr_free, self.ihead
  ptr_free, self.zhead
  ptr_free, self.ipsf
  ptr_free, self.zpsf
end

Pro Cluster::CreateCLEANCat, gals=gals
  if n_elements(gals) eq 0 then gals=indgen(n_elements(*self.galfitcat))
  if ~ptr_valid(self.CLEANcat) then begin
     CLEANcat = {CLEANcat, $
                 clusterid: '', $
                 galid: 0L, $
                 cleanap:  dblarr(5), $
                 cleanfac: dblarr(5) }
     CLEANcat = replicate(CLEANcat, n_elements(*self.sexcat))
     CLEANcat.clusterid = (*self.sexcat).clusterid
     CLEANcat.galid = (*self.sexcat).galid
     self.CLEANcat = ptr_new(CLEANcat)
  endif
  CLEANcat = *self.CLEANcat
  for igal=0, n_elements(gals)-1 do begin
     print, igal+1, ' of ', n_elements(gals)
     s1 = (self->summary())[gals[igal]]
     max_radius = 1.5*s1.r_e > 20.
     
     p = round([ s1.xwin_image-1, $ ;;sex coords to IDL coords
                 s1.ywin_image-1, $
                 s1.xwin_image-1, $
                 s1.ywin_image-1 ]) + [-1,-1,1,1]*round((max_radius*1.2) + 20)
     size=size(*self.iimage,/dim)
     if (p[0] lt 0 or p[1] lt 0 or p[2] ge size[0] or p[3] ge size[1]) then begin
        CLEANcat[gals[igal]].cleanap = dblarr(5)*!values.f_nan
        CLEANcat[gals[igal]].cleanfac = dblarr(5)*!values.f_nan
        continue
     endif
     iimage = self->ImageSection(p, band='i', /CLEAN)
     zimage = self->ImageSection(p, band='z', /CLEAN)
     apifluxes = apphot( iimage, $
                         (s1.xwin_image-0.5)-p[0], $
                         (s1.ywin_image-0.5)-p[1], $
                         [3.,5.,10.,15.,20.] )
     apzfluxes = apphot( zimage, $
                         (s1.xwin_image-0.5)-p[0], $
                         (s1.ywin_image-0.5)-p[1], $
                         [3.,5.,10.,15.,20.] )
     facifluxes = apphot( iimage, $
                          (s1.xwin_image-0.5)-p[0], $
                          (s1.ywin_image-0.5)-p[1], $
                          [0.5, 0.75, 1.0, 1.25, 1.5]*s1.r_e )
     faczfluxes = apphot( zimage, $
                          (s1.xwin_image-0.5)-p[0], $
                          (s1.ywin_image-0.5)-p[1], $
                          [0.5, 0.75, 1.0, 1.25, 1.5]*s1.r_e )
     isky = s1.isky
     zsky = s1.zsky
     if isky eq -1 then isky=s1.ibackground
     if zsky eq -1 then zsky=s1.zbackground
     apifluxes  -= isky*[3.,5.,10.,15.,20.]^2*!dpi
     apzfluxes  -= zsky*[3.,5.,10.,15.,20.]^2*!dpi
     facifluxes -= isky*([0.5, 0.75, 1.0, 1.25, 1.5]*s1.r_e)^2*!dpi
     faczfluxes -= zsky*([0.5, 0.75, 1.0, 1.25, 1.5]*s1.r_e)^2*!dpi
     CLEANcat[gals[igal]].cleanap = -2.5*alog10(apifluxes)+25.678 - (-2.5*alog10(apzfluxes)+24.867)
     CLEANcat[gals[igal]].cleanfac = -2.5*alog10(facifluxes)+25.678 - (-2.5*alog10(faczfluxes)+24.867)
  endfor
  ptr_free, self.CLEANcat
  self.CLEANcat = ptr_new(CLEANcat)
end

Pro Cluster__Define, struct
  struct = { Cluster, $
             Inherits JEMobject, $
             clustername:   '', $
             clusterid:     '', $
             ifilename:     '', $
             zfilename:     '', $
             zcluster:      0.d, $
             veldisp:       0.d, $
             ra:            0.d, $
             dec:           0.d, $
             ebv:           0.d, $
             sexcat:        ptr_new(), $
             galapagoscat:  ptr_new(), $
             skycat:        ptr_new(), $
             galfitcat:     ptr_new(), $
             morphcat:      ptr_new(), $
             LCcat:         ptr_new(), $
             CLEANcat:      ptr_new(), $
             speccat:       ptr_new(), $
             spectra:       ptr_new(), $
             iimage:        ptr_new(), $
             zimage:        ptr_new(), $
             ierrim:        ptr_new(), $
             zerrim:        ptr_new(), $
             ihead:         ptr_new(), $
             zhead:         ptr_new(), $
             iimepoch:      0L, $
             zimepoch:      0L, $
             ierepoch:      0L, $
             zerepoch:      0L, $
             maxepoch:      0L, $
             ipsf:          ptr_new(), $
             zpsf:          ptr_new(), $
             psffactor:     0L }
end
