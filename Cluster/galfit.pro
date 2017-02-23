Function Galfit, $
   iimage_ptr, $
   zimage_ptr, $
   ierrim_ptr, $
   zerrim_ptr, $
   galapagoscat, $
   clusterid, $
   ipsf, $
   zpsf, $
   gals=gals, $
   outdir=outdir, $
   oversample=oversample, $
   noerr=noerr, $
   iheader=iheader, $
   zheader=zheader

  if n_elements(oversample) eq 0 then oversample=1

  galfitcat = { galfitcat, $
                clusterid:'', $
                galid:0L, $
                ixpos:-1.d, $
                iypos:-1.d, $
                imag:-1.d, $
                iR_e:-1.d, $
                iboa:-1.d, $
                ipa:-1.d, $
                in:-1.d, $
                errixpos:-1.d, $
                erriypos:-1.d, $
                errimag:-1.d, $
                erriR_e:-1.d, $
                erriboa:-1.d, $
                erripa:-1.d, $
                errin:-1.d, $
                ichi2:-1.d, $
                idof:-1, $
                ichi2perdof:-1.d, $
                isky:-1.d, $
                isky_sigma:-1.d, $
                zxpos:-1.d, $
                zypos:-1.d, $
                zmag:-1.d, $
                zR_e:-1.d, $
                zboa:-1.d, $
                zpa:-1.d, $
                zn:-1.d, $
                errzxpos:-1.d, $
                errzypos:-1.d, $
                errzmag:-1.d, $
                errzR_e:-1.d, $
                errzboa:-1.d, $
                errzpa:-1.d, $
                errzn:-1.d, $
                zchi2:-1.d, $
                zdof:-1, $
                zchi2perdof:-1.d, $
                zsky:-1.d, $
                zsky_sigma:-1.d }

  galfitcat = replicate( galfitcat, n_elements(galapagoscat) )
  galfitcat.clusterid = galapagoscat.clusterid
  galfitcat.galid = galapagoscat.galid

  if n_elements(gals) eq 0 then gals = indgen( n_elements(galapagoscat) )
  

  ;;;;;;;;;;;;;;;
  ;; Loop through galaxies and fit

  for igal=0L, n_elements( gals ) - 1 do begin
     gal = galapagoscat[gals[igal]]

     ;;;;;;;;;;;;;;;;
     ;; BEGIN F775W
     psf_ptr = ptr_new(ipsf)
     gal_isophot_sky, iimage_ptr, ierrim_ptr, gal, galapagoscat, sky_value=sky_value, sky_sigma=sky_sigma
     galfitcat[gals[igal]].isky=sky_value
     galfitcat[gals[igal]].isky_sigma=sky_sigma
     galfit = MakeGalfit( iimage_ptr, ierrim_ptr, psf_ptr, $
                          gal, galapagoscat, sky_value, $
                          clusterid, 'i', $
                          outdir=outdir, oversample=oversample, $
                          noerr=noerr, iheader=iheader )
     ptr_free, psf_ptr
     galfit->Execute
     results = galfit->Results()
     if size(results, /tname) ne 'INT' then begin
        results2 = { iresults, $
                     ixpos:results.xpos + gal.xmin, $
                     iypos:results.ypos + gal.ymin, $
                     imag:results.mag, $
                     iR_e:results.R_e, $
                     iboa:results.boa, $
                     ipa:results.pa, $
                     in:results.n, $
                     errixpos:results.errxpos, $
                     erriypos:results.errypos, $
                     errimag:results.errmag, $
                     erriR_e:results.errR_e, $
                     erriboa:results.errboa, $
                     erripa:results.errpa, $
                     errin:results.errn, $
                     ichi2:results.chi2, $
                     idof:results.dof, $
                     ichi2perdof:results.chi2perdof }
        dest = galfitcat[gals[igal]]
        struct_assign, results2, dest, /nozero
        galfitcat[gals[igal]] = dest
     endif

     obj_destroy, galfit
     ;;;;;;;;;;;;;;;;
     ;; END F775W


     ;;;;;;;;;;;;;;;;
     ;; BEGIN F850LP
     psf_ptr = ptr_new(zpsf)
     gal_isophot_sky, zimage_ptr, zerrim_ptr, gal, galapagoscat, sky_value=sky_value, sky_sigma=sky_sigma
     galfitcat[gals[igal]].zsky=sky_value
     galfitcat[gals[igal]].zsky_sigma=sky_sigma
     galfit = MakeGalfit( zimage_ptr, zerrim_ptr, psf_ptr, $
                          gal, galapagoscat, sky_value, $
                          clusterid, 'z', $
                          outdir=outdir, oversample=oversample, $
                          noerr=noerr, zheader=zheader )
     ptr_free, psf_ptr
     galfit->Execute
     results = galfit->Results()
     if size(results, /tname) ne 'INT' then begin
        results2 = { zresults, $
                     zxpos:results.xpos+gal.xmin, $
                     zypos:results.ypos+gal.ymin, $
                     zmag:results.mag, $
                     zR_e:results.R_e, $
                     zboa:results.boa, $
                     zpa:results.pa, $
                     zn:results.n, $
                     errzxpos:results.errxpos, $
                     errzypos:results.errypos, $
                     errzmag:results.errmag, $
                     errzR_e:results.errR_e, $
                     errzboa:results.errboa, $
                     errzpa:results.errpa, $
                     errzn:results.errn, $
                     zchi2:results.chi2, $
                     zdof:results.dof, $
                     zchi2perdof:results.chi2perdof }
        dest = galfitcat[gals[igal]]
        struct_assign, results2, dest, /nozero
        galfitcat[gals[igal]] = dest
     endif

     obj_destroy, galfit
     ;;;;;;;;;;;;;;;;
     ;; END F850LP

  endfor

  return, galfitcat
end
