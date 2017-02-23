Function MakeGalfit1, $
   image_ptr, $
   errim_ptr, $
   psf_ptr, $
   galapagoscat, $
   sexcat, $
   skycat, $
   igal, $
   band, $
   iband, $
   zeropoint, $
   realsky, $
   outdir=outdir, $
   header=header, $
   noerr=noerr, $
   take2=take2, $
   snx=snx, $
   sny=sny, $
   snmag=snmag

  if ~finite(realsky) then begin
     realsky=0.
     if band eq 'z' then realsky = 0.0268
     if band eq 'i' then realsky = 0.0504
  endif

  psf=*psf_ptr
  ;; grab image data
  nx = (size( (*image_ptr), /dimens ))[0]
  ny = (size( (*image_ptr), /dimens ))[1]

  gal1=galapagoscat[igal]
  sex1=sexcat[igal,iband]
  sky1=skycat[igal,iband]
  xovermax = gal1.xmax-(nx-1)
  yovermax = gal1.ymax-(ny-1)
  if xovermax gt 0 then begin
     gal1.xmax -= xovermax
     gal1.xmin -= xovermax
  endif
  if yovermax gt 0 then begin
     gal1.ymax -= yovermax
     gal1.ymin -= yovermax
  endif
  if gal1.xmin lt 0 then return, 0
  if gal1.ymin lt 0 then return, 0

  mag=sex1.mag_auto
  re=(0.162d)*sex1.flux_radius^(1.87d)
  sky=sky1.sky_value
  if ~finite(sky) then return, 0
  xpos=sex1.xwin_image-gal1.xmin
  ypos=sex1.ywin_image-gal1.ymin
  n=1.5
  boa=sex1.b_image/sex1.a_image
  pa=sex1.theta_image-90.

  params = { psffactor:   1, $
             fit_xmin:    1, $
             fit_ymin:    1, $
             fit_xmax:    gal1.xmax-gal1.xmin+1, $
             fit_ymax:    gal1.ymax-gal1.ymin+1, $
             conv_x:      45, $
             conv_y:      45, $
             display:     'regular', $
             interactive: 0, $
             zeropoint:   zeropoint, $
             platescalex: 0.05, $
             platescaley: 0.05, $
             options:     0 }
  
  constraint1 = obj_new( 'galfitabsconstraint', 1, 're', 0.18, 300. )
  constraint2 = obj_new( 'galfitrelconstraint', 1, 'mag', -5., 5. )
  constraint3 = obj_new( 'galfitabsconstraint', 1, 'n', 1.0d, 4.0d )

  constraints = [constraint1, constraint2, constraint3]
  
  component1 = obj_new( 'galfitsersic', xpos=xpos, ypos=ypos, $
                        mag=mag, re=re, n=n, boa=boa, pa=pa, out=0)
  components = component1

  ;; add any simultaneously fit objects
  if (~keyword_set(take2) and ptr_valid(gal1.simfit)) $
     or (keyword_set(take2) and ptr_valid(gal1.simfit2)) then begin

     if keyword_set(take2) then simfit = *gal1.simfit2 $
     else simfit = *gal1.simfit

     for isim=0L, n_elements(simfit)-1 do begin
        sex2  = sexcat[simfit[isim],iband]
        sxpos = sex2.xwin_image-gal1.xmin
        sypos = sex2.ywin_image-gal1.ymin
        smag  = sex2.mag_auto
        sre   = (0.162d)*(sex2.flux_radius)^(1.87d)
        sn    = 1.5
        sboa  = sex2.b_image/sex2.a_image
        spa   = sex2.theta_image-90.
        component = obj_new( 'galfitsersic', xpos=sxpos, ypos=sypos, $
                             mag=smag, re=sre, n=sn, boa=sboa, pa=spa, out=0 )
        components = [components, component]
        constraint1 = obj_new( 'galfitabsconstraint', isim+2, 'n', 0.5d, 8.0d )
        constraint2 = obj_new( 'galfitabsconstraint', isim+2, 're', 0.18, 600 )
        constraints = [constraints, constraint1, constraint2]
     endfor
  endif

  ;; add supernova PSF if requested
  if n_elements(snx) ne 0 then begin
     snxpos = snx-gal1.xmin
     snypos = sny-gal1.ymin
     snmag  = snmag
     component = obj_new( 'galfitpsf', xpos=snxpos, ypos=snypos, mag=snmag, out=0 )
     components = [components, component]
  endif
  
  ;; add sky as an object
  if keyword_set(noerr) then $
     isophotsky = (sky+realsky)*sxpar(header, 'EXPTIME') $
  else $
     isophotsky = sky
  skycomponent = obj_new( 'galfitsky', sky=isophotsky, out=0)
  components = [components, skycomponent]
  
  ;; create mask for other interloping objects
  pnx = gal1.xmax-gal1.xmin+1
  pny = gal1.ymax-gal1.ymin+1
  mskim = dblarr( pnx, pny )
  if (~keyword_set(take2) and ptr_valid(gal1.maskfit)) $
     or (keyword_set(take2) and ptr_valid(gal1.maskfit2)) then begin

     if keyword_set(take2) then maskfit = *gal1.maskfit2 $
     else maskfit = *gal1.maskfit

     for imask=0L, n_elements(maskfit)-1 do begin
        sex3 = sexcat[maskfit[imask],0]
        gal3 = galapagoscat[maskfit[imask]]
        dist_ellipse, ellipse, [pnx, pny], $
                      sex3.xwin_image-gal1.xmin-1, sex3.ywin_image-gal1.ymin-1, $
                      sex3.a_image/sex3.b_image, sex3.theta_image-90.
        w = where( ellipse lt gal3.majaxis )
        if w[0] ne -1 then mskim[w] = 1
     endfor
  endif

  prefix = gal1.clusterid+string( gal1.galid, format='(I05)' )+band
  image = (*image_ptr)[gal1.xmin:gal1.xmax, gal1.ymin:gal1.ymax]
  errim = (*errim_ptr)[gal1.xmin:gal1.xmax, gal1.ymin:gal1.ymax]
  mskim = mskim or ~finite(errim)
  mkhdr, hdr, image ;;need to redo header to get right image size...
  sxaddpar, hdr, 'EXPTIME', sxpar(header, 'EXPTIME')
  sxaddpar, hdr, 'GAIN', sxpar(header, 'GAIN')
  sxaddpar, hdr, 'RDNOISE', sxpar(header, 'RDNOISE')
  sxaddpar, hdr, 'NCOMBINE', sxpar(header, 'NCOMBINE')
  if keyword_set(noerr) then begin
     image += (realsky+sky)
     image *= sxpar(header, 'EXPTIME')
     galfit = obj_new( 'galfit', params=params, components=components, $
                       image=image, mskim=mskim, psf=psf, constraints=constraints, $
                       outdir=outdir, prefix=prefix, header=hdr )
  endif else begin
     galfit = obj_new( 'galfit', params=params, components=components, $
                       image=image, errim=errim, mskim=mskim, $
                       psf=psf, constraints=constraints, $
                       outdir=outdir, prefix=prefix )
  endelse
  return, galfit
end
