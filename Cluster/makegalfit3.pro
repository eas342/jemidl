Function MakeGalfit3, sciim, $
                      errim, $
                      psf, $
                      sexcat, $   ;;single band
                      galapagoscat, $
                      skycat, $   ;;single band
                      igal, $
                      hdr, $
                      outdir=outdir, $
                      sn=sn, $
                      take2=take2, $
                      noerr=noerr

  sz = size( sciim, /dim )
  gal1 = galapagoscat[igal]
  sex1 = sexcat[igal]
  sky1 = skycat[igal]

  params = { psffactor:   1, $
             fit_xmin:    1, $
             fit_ymin:    1, $
             fit_xmax:    gal1.xmax-gal1.xmin+1, $
             fit_ymax:    gal1.ymax-gal1.ymin+1, $
             conv_x:      45, $
             conv_y:      45, $
             display:     'regular', $
             interactive: 0, $
             zeropoint:   hdr.zeropoint, $
             platescalex: 0.05, $
             platescaley: 0.05, $
             options:     0 }
  
  ;;target object
  component1 = obj_new( 'galfitsersic', $
                        xpos=sex1.xwin_image-gal1.xmin, $
                        ypos=sex1.ywin_image-gal1.ymin, $
                        mag=sex1.mag_auto, $
                        re=(0.162d)*sex1.flux_radius^(1.87d), $
                        n=1.5, $
                        boa=sex1.b_image/sex1.a_image, $
                        pa=sex1.theta_image-90., $
                        out=0 )
  constraint1 = obj_new( 'galfitabsconstraint', 1, 're', 0.18, 300.0 )
  constraint2 = obj_new( 'galfitrelconstraint', 1, 'mag', -5.0, 5.0 )
  constraint3 = obj_new( 'galfitabsconstraint', 1, 'n', 1.0d, 4.0d )
  constraints=[constraint1, constraint2, constraint3]
  components=component1

  ;; add any simultaneously fit objects and create mask
  if (~keyword_set(take2) and ptr_valid(gal1.simfit)) $
     or (keyword_set(take2) and ptr_valid(gal1.simfit2)) then begin

     if keyword_set(take2) then begin
        simfit = *gal1.simfit2
     endif else begin
        simfit = *gal1.simfit
     endelse

     for isim=0L, n_elements(simfit)-1 do begin
        sex2 = sexcat[simfit[isim]]
        component = obj_new( 'galfitsersic', $
                             xpos=sex2.xwin_image-gal1.xmin, $
                             ypos=sex2.ywin_image-gal1.ymin, $
                             mag=sex2.mag_auto, $
                             re=0.162d*sex2.flux_radius^(1.87d), $
                             n=1.5, $
                             boa=sex2.b_image/sex2.a_image, $
                             pa=sex2.theta_image-90.0, $
                             out=0 )
        components = [components, component]
        constraint1 = obj_new( 'galfitabsconstraint', isim+2, 'n', 0.5d, 8.0d )
        constraint2 = obj_new( 'galfitabsconstraint', isim+2, 're', 0.18, 600 )
        constraints = [constraints, constraint1, constraint2]
     endfor
  endif

  ;;make mask
  mskim = bytarr( sz[0], sz[1] )
  if (~keyword_set(take2) and ptr_valid(gal1.maskfit)) $
     or (keyword_set(take2) and ptr_valid(gal1.maskfit2)) then begin

     if keyword_set(take2) then begin
        maskfit = *gal1.maskfit2
     endif else begin
        maskfit = *gal1.maskfit
     endelse

     for imask=0L, n_elements(maskfit)-1 do begin
        sex3 = sexcat[maskfit[imask]]
        gal3 = galapagoscat[maskfit[imask]]
        dist_ellipse, ellipse, sz, $
                      sex3.xwin_image-gal1.xmin-1, sex3.ywin_image-gal1.ymin-1, $
                      sex3.a_image/sex3.b_image, sex3.theta_image-90.
        w = where( ellipse lt gal3.majaxis )
        if w[0] ne -1 then mskim[w] = 1
     endfor
  endif

  ;;add supernova PSF if requested
  if n_elements(sn) ne 0 then begin
     component = obj_new( 'galfitpsf', $
                          xpos=sn.x-gal1.xmin, $
                          ypos=sn.y-gal1.ymin, $
                          mag=sn.mag, $
                          out=0 )
     components = [components, component]
  endif
  
  ;;add sky as an object
  if ~finite(sky1.sky_value) then return, -1
  if keyword_set(noerr) then $
     skyval = (sky1.sky_value+sky1.skyim_value)*hdr.exptime $
  else $
     skyval = sky1.sky_value
  skycomponent = obj_new( 'galfitsky', sky=skyval, out=0 )
  components = [components, skycomponent]

  ;;create header
  prefix = gal1.clusterid+string( gal1.galid, format='(I05)' )+hdr.band
  mskim = mskim or ~(finite(errim))
  mkhdr, hdr0, sciim
  sxaddpar, hdr0, 'EXPTIME', hdr.exptime
  sxaddpar, hdr0, 'GAIN', hdr.gain
  sxaddpar, hdr0, 'RDNOISE', hdr.rdnoise
  sxaddpar, hdr0, 'NCOMBINE', hdr.ncombine

  ;;put it all together
  if keyword_set(noerr) then begin
     sciim += (sky1.sky_value+sky1.skyim_value)
     sciim *= hdr.exptime
     galfit = obj_new( 'galfit', params=params, components=components, $
                       image=sciim, mskim=mskim, $
                       psf=psf, constraints=constraints, $
                       outdir=outdir, prefix=prefix, header=hdr0 )
  endif else begin
     galfit = obj_new( 'galfit', params=params, components=components, $
                       image=sciim, errim=errim, mskim=mskim, $
                       psf=psf, constraints=constraints, $
                       outdir=outdir, prefix=prefix )
  endelse
  return, galfit
end
