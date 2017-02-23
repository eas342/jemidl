Function hdr2helio, hdr
  telescope = strcompress(strmid(sxpar(hdr[*, 0], 'TELESCOP'), 0, 3) $
                          , /rem)
  IF NOT KEYWORD_SET(TELESCOPE) THEN TELESCOPE = $
     strcompress(strmid(sxpar(hdr[*, 0], 'TELID'), 0, 3), /rem)
  if telescope eq '' and strcompress(sxpar(hdr[*,0], 'INSTRUME'),/rem) eq 'OSIRIS' then telescope='GTC'

  case telescope OF
     'Kec': begin  ;; Keck/LRIS
        mjd = double(sxpar(hdr, 'MJD-OBS'))  + 2400000.5D
        equinox = double(sxpar(hdr, 'EQUINOX'))
        ra = sxpar(hdr, 'RA')
        dec = sxpar(hdr, 'DEC')
        x_radec, ra, dec, radeg, decdeg
     end
     'GTC': begin ;; GTC
        mjd = sxpar(hdr, 'MJD-OBS')+2400000.5D
        equinox = 2000.0
        obs = 'lapalma'
        ra = sxpar(hdr, 'RA')
        dec = sxpar(hdr, 'DEC')
        x_radec, ra, dec, radeg, decdeg
     end
  endcase
  if not keyword_set(UNDO) then SIGN = -1. else SIGN = 1.
  helio = (SIGN)*x_keckhelio(radeg, decdeg, equinox, jd = mjd, OBS = OBS)
  hel_corr = sqrt( (1.d + helio/299792.458d) / $
                   (1.d - helio/299792.458d) )
  return, hel_corr ;; km/s
end

Pro coaddscience, $
   dir, $
   sensfuncfile, $
   waverange, $
   objidfile=objidfile, $
   tellfile=tellfile, $
   check=check, $
   CR=CR

  dir=directoryify(dir)
  if keyword_set(CR) then begin
     outdir='FinalCR/'
  endif else begin
     outdir='Final/'
  endelse
  if strpos(sensfuncfile,'bsens') ne -1 then $
     searchstr = 'sci-*.fits.gz' $
  else $
     searchstr = 'sci-*.newwave.fits.gz'

  f=file_search(dir+searchstr, count=nfiles)
  targs = strarr(nfiles)
  for i=0, nfiles-1 do begin
     hdr=headfits(f[i], /silent)
     if strmatch(strtrim(sxpar(hdr,'INSTRUME'),2),'OSIRIS') then begin
        targs[i]=lowzhost_nametranslate(strtrim(sxpar(hdr, 'OBJECT'), 2))
     endif else begin
        targs[i]=lowzhost_nametranslate(strtrim(sxpar(hdr, 'TARGNAME'), 2))
     endelse
     print, 'loading: '+f[i]+' '+targs[i]
  endfor
  print
  uniqtargs = targs[uniq(targs,sort(targs))]
  for i=0, n_elements(uniqtargs)-1 do begin
     name = uniqtargs[i]
     print, name
     if ~(file_info(outdir)).directory then begin
        file_mkdir, outdir
     endif
     nfil = outdir + name + '_N.fits' ;; coadded, not fluxed
     ffil = outdir + name + '_F.fits' ;; fluxed (optionally telluric-corrected)
     rfil = outdir + name + '_R.fits' ;; resolution
     ntfil = outdir + name + '_NT.fits' ;; not telluric-corrected
     chekfile = outdir + name + '_check.fits'
     bpixfile = 'badpix/' + name + '.pix'
     if ~(file_info(bpixfile)).exists or ~keyword_set(CR) then bpixfile=''

     if (file_info(nfil)).exists then continue ;; don't overwrite old results?

     delvarx, infiles, objid
     for j=0, n_elements(targs)-1 do begin
        if name ne targs[j] then continue
        print, ' ', f[j]
        if n_elements(infiles) eq 0 then begin
           infiles = f[j]
        endif else begin
           infiles = [infiles, f[j]]
        endelse
     endfor

     ;; find appropriate objects
     objid = fltarr(n_elements(infiles))
     ;; first check for existence of objidfile
     if n_elements(objidfile) ne 0 then begin
        readcol, objidfile, targname, objidstr, f='A,A', delim=' ', /silent
        w=where(name eq targname)
        if w[0] eq -1 then message, "Couldn't find targname: "+name+" in objidfile"
        objid = long(strsplit(objidstr[w[0]], ',', /extract))
     endif else begin
        ;; if not objidfile specified, pick brightest object on slit for each exposure
        for j=0, n_elements(infiles)-1 do begin
           a=mrdfits(infiles[j],5)
           nspec=n_elements(a)
           bright = fltarr(nspec)
           for k=0, nspec-1 do begin
              bright[k] = median(a[k].flux_opt)
           endfor
           junk = max(bright,m)
           objid[j] = m+1
        endfor
     endelse

     jm_long_coadd, infiles, objid, outfil=nfil, check=check, chekfile=chekfile, $
                    wave=waven, flux=fluxn, ivar=ivarn, /box, badpixfile=bpixfile
     long_fluxcal, nfil, sensfuncfile=sensfuncfile, outfil=ffil, wave=wave, flux=flux, sig=sig
     if n_elements(tellfile) ne 0 then begin
        tell=mrdfits(tellfile, 0, shdr, /silent)
        twave=mrdfits(tellfile, 1, /silent)
        std_airmass = sxpar(shdr, 'AIRMASS')
        fhdr = headfits(ffil)
        obj_airmass = sxpar(fhdr, 'AIRMASS')
        influx = x_readspec(ffil, wav=wave, sig=sig, inflg=2)
        ;; need to compensate for helio-shifts
        std_hcorr = hdr2helio(shdr)
        obj_hcorr = hdr2helio(fhdr)
        tell = interpol(tell, twave/std_hcorr, wave/obj_hcorr, /spline)
        tflux = influx / (tell^((obj_airmass/std_airmass)^0.6)) ;; 0.6 from Buton2012
        tsig = sig / (tell^((obj_airmass/std_airmass)^0.6))
        file_move, ffil, ntfil, /overwrite
        mwrfits, tflux, ffil, fhdr, /create
        mwrfits, tsig, ffil
        mwrfits, wave, ffil
     endif
     influx = x_readspec(ffil, wav=wave, sig=sig, inflg=2)
     long_resvec, infiles, objid, wave, resvec=resvec, outfil=rfil
     inivar = (sig gt 0.0)/(sig^2 + (sig eq 0.0))
     inloglam = alog10(wave)
     imask = (wave le waverange[0]) or (wave ge waverange[1])
     ibad = where(imask)
     inivar[ibad] = 0.0d
     gal_struct = jm_long_galstats(influx, inivar, resvec, inloglam, $
                                   inter=check, /telluric, /splinemask, $
                                   vdispmodel=vdispmodel, hardcopy=outdir+name)
     if keyword_set(check) then x_specplot, ffil, inflg=2, xsiz=1200, ysiz=600, $
                                            block=check, zin=gal_struct[0].z, $
                                            /gal, ytwo=(imask eq 0)*vdispmodel, $
                                            two_wave=wave
     mwrfits, gal_struct, outdir+name+'_galstruct.fits'
     lickcheck, chekfile
  endfor
end
