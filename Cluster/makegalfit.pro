Function MakeGalfit, $
   image_ptr, $
   errim_ptr, $
   psf_ptr, $
   gal, $
   galapagoscat, $
   band, $
   outdir=outdir, $
   header=header, $
   noerr=noerr, $
   take2=take2

  psf=*psf_ptr
  ;; grab image data
  nx = (size( (*image_ptr), /dimens ))[0]
  ny = (size( (*image_ptr), /dimens ))[1]

  xovermax = gal.xmax-(nx-1)
  yovermax = gal.ymax-(ny-1)
  if xovermax gt 0 then begin
     gal.xmax -= xovermax
     gal.xmin -= xovermax
  endif
  if yovermax gt 0 then begin
     gal.ymax -= yovermax
     gal.ymin -= yovermax
  endif

  case band of 
     'i' : begin
        mag_guess=gal.imag_guess
        zp=25.678
        r_e_guess=gal.ir_e_guess
        sky=gal.isky
        realsky = 0.0504
     end
     'z' : begin
        mag_guess=gal.zmag_guess
        zp=24.867
        r_e_guess=gal.zr_e_guess
        sky=gal.zsky
        realsky = 0.0268
     end
  endcase

  params = { psffactor:   1, $
             fit_xmin:    1, $
             fit_ymin:    1, $
             fit_xmax:    gal.xmax-gal.xmin+1, $
             fit_ymax:    gal.ymax-gal.ymin+1, $
             conv_x:      45, $
             conv_y:      45, $
             display:     'regular', $
             interactive: 0, $
             zeropoint:   zp, $
             platescalex: 0.05, $
             platescaley: 0.05, $
             options:     0 }
  
  constraint1 = obj_new( 'galfitabsconstraint', 1, 're', 0.18, 300. )
  constraint2 = obj_new( 'galfitrelconstraint', 1, 'mag', -5., 5. )
  constraint3 = obj_new( 'galfitabsconstraint', 1, 'n', 1.0d, 4.0d )

  constraints = [constraint1, constraint2, constraint3]
  
  component1 = obj_new( 'galfitsersic', xpos=gal.x_guess-gal.xmin, ypos=gal.y_guess-gal.ymin, $
                        mag=mag_guess, r_e=r_e_guess, n=gal.n_guess, boa=gal.boa_guess, $
                        pa=gal.pa_guess, out=0)
  components = component1

  ;; add any simultaneously fit objects
  if (~keyword_set(take2) and ptr_valid(gal.simfit)) or (keyword_set(take2) and ptr_valid(gal.simfit2)) then begin
     if keyword_set(take2) then simfit = *gal.simfit2 $
     else simfit = *gal.simfit

     for isim=0L, n_elements(simfit)-1 do begin
        simgal = galapagoscat[(simfit)[isim]]
        case band of 
           'i' : begin
              smag_guess=simgal.imag_guess
              sr_e_guess=simgal.ir_e_guess
           end
           'z' : begin
              smag_guess=simgal.zmag_guess
              sr_e_guess=simgal.zr_e_guess
           end
        endcase
        component = obj_new( 'galfitsersic', xpos=simgal.x_guess-gal.xmin, ypos=simgal.y_guess-gal.ymin, $
                             mag=smag_guess, r_e=sr_e_guess, n=simgal.n_guess, $
                             boa=simgal.boa_guess, pa=simgal.pa_guess, out=0 )
        components = [components, component]
     endfor
  endif
  
  ;; add sky as an object
  if keyword_set(noerr) then $
     isophotsky = (sky+realsky)*sxpar(header, 'EXPTIME') $
  else $
     isophotsky = sky
  skycomponent = obj_new( 'galfitsky', sky=isophotsky, out=0)
  components = [components, skycomponent]
  
  ;; create mask for other interloping objects
  pnx = gal.xmax-gal.xmin+1
  pny = gal.ymax-gal.ymin+1
  mskim = dblarr( pnx, pny )
  if (~keyword_set(take2) and ptr_valid(gal.maskfit)) or (keyword_set(take2) and ptr_valid(gal.maskfit2)) then begin
     if keyword_set(take2) then maskfit = *gal.maskfit2 $
     else maskfit = *gal.maskfit

     for imask=0L, n_elements(maskfit)-1 do begin
        maskgal = galapagoscat[(maskfit)[imask]]
        dist_ellipse, ellipse, [pnx, pny], $
                      maskgal.x_guess-gal.xmin-1, maskgal.y_guess-gal.ymin-1, $
                      1./maskgal.boa_guess, maskgal.pa_guess
        w = where( ellipse lt maskgal.majaxis )
        if w[0] ne -1 then mskim[where( ellipse lt maskgal.majaxis )] = 1
     endfor
  endif
  
  prefix = gal.clusterid+string( gal.galid, format='(I05)' )+band
  image = (*image_ptr)[gal.xmin:gal.xmax, gal.ymin:gal.ymax]
  errim = (*errim_ptr)[gal.xmin:gal.xmax, gal.ymin:gal.ymax]
  mkhdr, hdr, image
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
                       image=image, $
                       errim=errim, $
                       mskim=mskim, psf=psf, constraints=constraints, $
                       outdir=outdir, prefix=prefix )
  endelse
  return, galfit
end
