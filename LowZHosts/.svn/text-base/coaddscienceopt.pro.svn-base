Pro coaddscienceopt, dir, sensfuncfile, waverange, check=check
  dir=directoryify(dir)
  if strpos(sensfuncfile,'bsens') ne -1 then $ ;; blue or red?
     searchstr = 'sci-*.fits.gz' $
  else $
     searchstr = 'sci-*.newwave.fits.gz'

  f=file_search(dir+searchstr, count=nfiles)
  targs = strarr(nfiles)
  for i=0, nfiles-1 do begin
     hdr=headfits(f[i], /silent)
     targs[i]=lowzhost_nametranslate(strtrim(sxpar(hdr, 'TARGNAME'), 2))
  endfor
  uniqtargs = targs[uniq(targs,sort(targs))]
  for i=0, n_elements(uniqtargs)-1 do begin
     name = uniqtargs[i]
     print, name
     if ~(file_info('FinalOpt/')).directory then begin
        file_mkdir, 'FinalOpt'
     endif
     nfil = 'FinalOpt/' + name + '_N.fits'
     ffil = 'FinalOpt/' + name + '_F.fits'
     rfil = 'FinalOpt/' + name + '_R.fits'
     chekfile = 'FinalOpt/' + name + '_check.fits'

     if (file_info(nfil)).exists then continue

     delvarx, infiles, objid
     for j=0, n_elements(targs)-1 do begin
        if name ne targs[j] then continue
        print, ' ', f[j]
        if n_elements(infiles) eq 0 then begin
           infiles = f[j]
        endif else begin
           infiles = [infiles, f[j]]
        endelse
        objid = replicate(1, n_elements(infiles))
     endfor
     jm_long_coadd, infiles, objid, outfil=nfil, check=check, chekfile=chekfile, $
                    wave=waven, flux=fluxn, ivar=ivarn
     long_fluxcal, nfil, sensfuncfile=sensfuncfile, outfil=ffil, wave=wave, flux=flux, sig=sig
     influx = x_readspec(ffil, wav=wave, sig=sig, inflg=2)
     long_resvec, infiles, objid, wave, resvec=resvec, outfil=rfil
     inivar = (sig gt 0.0)/(sig^2 + (sig eq 0.0))
     inloglam = alog10(wave)
     imask = (wave le waverange[0]) or (wave ge waverange[1])
     ibad = where(imask)
     inivar[ibad] = 0.0d
     gal_struct = jm_long_galstats(influx, inivar, resvec, inloglam, $
                                   inter=check, /telluric, /splinemask, $
                                   vdispmodel=vdispmodel, hardcopy='FinalOpt/'+name)
     if keyword_set(check) then x_specplot, ffil, inflg=2, xsiz=1200, ysiz=600, $
                                            block=check, zin=gal_struct[0].z, $
                                            /gal, ytwo=(imask eq 0)*vdispmodel, $
                                            two_wave=wave
     mwrfits, gal_struct, 'FinalOpt/'+name+'_galstruct.fits'
     lickcheck, chekfile
  endfor
end
