Pro indiv_index_check, chekfile, setup
  set_plot, 'PS'
  loadct, 0
  skyblue   = fsc_color('skyblue',   255)
  papaya    = fsc_color('papaya',    254)
  darkgreen = fsc_color('darkgreen', 253)
  red       = fsc_color('red',       252)

  ;; Load Lick index wavelength definitions
  lickfile=getenv('EZ_AGES_DIR')+'lickindexlist.txt'
  readcol, lickfile, bp0, bp1, bc0, bc1, rc0, rc1, lname, $
           format='F,F,F,F,F,F,X,X,A', /silent

  ;; Adding some emission lines by hand
  ;; bp = bandpass
  ;; bc = blue continuum
  ;; rc = red continuum
  bp0   = [bp0,   [  3719.65,    4942.50,    4991.25]]
  bp1   = [bp1,   [  3737.35,    4980.50,    5029.25]]
  bc0   = [bc0,   [  3660.00,    4876.70,    4876.70]]
  bc1   = [bc1,   [  3699.00,    4933.60,    4933.60]]
  rc0   = [rc0,   [  3756.35,    5051.50,    5051.50]]
  rc1   = [rc1,   [  3789.30,    5118.60,    5118.60]]
  lname = [lname, ['OII3727', 'OIII4959', 'OIII5007']]

  ;; And some more
  bp0   = [bp0,   [  6542.65,  6556.85,   6577.25,   6712.65,   6726.40]]
  bp1   = [bp1,   [  6556.85,  6572.25,   6592.65,   6724.50,   6739.00]]
  bc0   = [bc0,   [  6515.00,  6515.00,   6515.00,   6685.00,   6685.00]]
  bc1   = [bc1,   [  6538.00,  6538.00,   6538.00,   6706.00,   6706.00]]
  rc0   = [rc0,   [  6600.00,  6600.00,   6600.00,   6744.00,   6744.00]]
  rc1   = [rc1,   [  6625.00,  6625.00,   6625.00,   6764.00,   6764.00]]
  lname = [lname, ['NII6548', 'Halpha', 'NII6584', 'SII6717', 'SII6731']]

  ;; Use nice latex names for plot titles
  lnametex = repstr(lname, 'alpha', '\alpha')
  lnametex = repstr(lnametex, 'beta', '\beta')
  lnametex = repstr(lnametex, 'gamma', '\gamma')
  lnametex = repstr(lnametex, 'delta', '\delta')
  lnametex = repstr(lnametex, 'OII3', '[OII]3') ;; OII and OIII interfere without the 3
  lnametex = repstr(lnametex, 'OIII', '[OIII]')
  lnametex = repstr(lnametex, 'NII', '[NII]')
  lnametex = repstr(lnametex, 'SII', '[SII]')

  check = mrdfits(chekfile, 1, /silent)
  sz = size(check.indiv, /dim)
  nobj = sz[0]
  npix = sz[1]
  SNname = repstr(chekfile, '_check.fits', '')
  galstruct = mrdfits(SNname+'_galstruct.fits', 1, /silent)
  fluxfile = SNname+'_F.fits'

  h = headfits(fluxfile, /silent)
  flux = x_readspec(fluxfile, sig=sig, inflg=2)
  S2N = flux/sig
  wave = check.wave/(1.0+galstruct[0].z)

  ;; Open pixel flat
  dir = strsplit('./'+chekfile, '/', /extract)
  dir = strjoin(dir[0:n_elements(dir)-2], '/', /single)+'/'
  pixflat = mrdfits((file_search(dir+'../pixflat*fits'))[0], /silent)

  ;; Blue data is flipped
  if member(setup, ['B460','B560']) then $
     pixflat = reverse(pixflat, 2)

  ;; Use long SN names instead of telescope target names for output files
  SNoutname = lowzhost_nametranslate(SNname)

  ;; Open Science data
  objstruct = ptrarr(nobj)
  scistruct = ptrarr(nobj)
  objidfile = '../objid'
  if (file_info(objidfile)).exists then begin
     readcol, objidfile, targname, objidstr, f='A,A', delim=' ', /silent
     objid = fix(strsplit(objidstr, ',', /extract))
  endif
  for iobj=0, nobj-1 do begin
     scifile = dir+'../Science_box/'+sxpar(h, 'INFILE'+strn(iobj,f='(I02)'))
     obj = mrdfits(scifile, 5, /silent)
     ntrace = n_elements(obj)
     signal = fltarr(ntrace)
     if n_elements(objid) eq 0 then begin
        for itrace=0, ntrace-1 do begin
           signal[itrace] = median(obj[itrace].flux_box)
        endfor
        junk = max(signal,m)
     endif else begin
        w=where(targname eq SNoutname)
        if w[0] eq -1 then message, 'SN not found in objidfile'
        m = objid[iobj]-1
     endelse
     objstruct[iobj] = ptr_new(obj[m])
     sci = mrdfits(repstr(scifile, '.newwave', ''), 0, /silent)
     sky = mrdfits(repstr(scifile, '.newwave', ''), 2, /silent)
     if member(setup, ['B460','B560']) then $
        scistruct[iobj] = ptr_new(reverse(sci-sky,2)) $
     else $
        scistruct[iobj] = ptr_new(sci-sky)
  endfor

  if ~(file_info('LickEPS/')).exists then file_mkdir, 'LickEPS'
  for i=0, n_elements(lname)-1 do begin ;; loop over Lick indices
     xrange = [bc0[i], rc1[i]]
     w=where(wave ge xrange[0] and wave le xrange[1])
     if w[0] eq -1 then continue ;; no data in wavelength range

     medS2N = median(S2N[w])
     w1 = where(~check.mask $
                and check.indiv ne 0 $
                and transpose(rebin(wave,npix,nobj)) ge xrange[0] $
                and transpose(rebin(wave,npix,nobj)) le xrange[1])
     yrange = percentile(check.indiv[w1],[0.02,0.98])
     xrange += 0.1*[-1,1]*(xrange[1]-xrange[0])
     yrange += 0.3*[-1,1]*(yrange[1]-yrange[0])

     ;; First plot is coadd with all indivs overplotted
     ps_start, filename='LickEPS/'+SNoutname+'.'+setup+'.'+lname[i]+'.coadd.eps'
     plot, [0], /nodata, xrange=xrange, yrange=yrange, $
           xstyle=5, ystyle=1, position=[0.12, 0.10, 0.97, 0.90]
;     plot, [0], /nodata, xrange=xrange, yrange=yrange, $
;           xstyle=5, ystyle=1, xmargin=[7,2], ymargin=[3,4]
     xyouts, 0.5, 0.95, SNoutname, align=0.5, /normal
     if finite(medS2N) then $
        xyouts, 0.06, 0.95, textoidl(lnametex[i])+' S/N ~ '+strn(medS2N,F='(I)'), $
                align=0, /normal $
     else $
        xyouts, 0.06, 0.95, textoidl(lnametex[i]), align=0, /normal
     xyouts, 0.97, 0.95, setup+'/'+'COADD', $
             align=1, /normal, charsize=1.0
     polyfill, [bc0[i], bc0[i], bc1[i], bc1[i], bc0[i]], $
               [yrange[0], yrange[1], yrange[1], yrange[0], yrange[0]], $
               color=skyblue
     polyfill, [rc0[i], rc0[i], rc1[i], rc1[i], rc0[i]], $
               [yrange[0], yrange[1], yrange[1], yrange[0], yrange[0]], $
               color=skyblue
     polyfill, [bp0[i], bp0[i], bp1[i], bp1[i], bp0[i]], $
               [yrange[0], yrange[1], yrange[1], yrange[0], yrange[0]], $
               color=papaya
     axis, xaxis=1, xrange=xrange*(1.0+galstruct[0].z), xstyle=1 ;; obs wavelength
     axis, xaxis=0, xrange=xrange, xstyle=1 ;; rest wavelength
     oplot, wave, check.coadd, thick=6, ps=10
     w1 = where(check.indiv eq 0.0) ;; don't plot fluxes of 0.0
     if w1[0] ne -1 then check.indiv[w1] = !values.f_nan
     for iobj=0, nobj-1 do begin
        oplot, wave, check.indiv[iobj,*], thick=2, ps=10, color=darkgreen
        wbad = where(check.mask[iobj,*], nbad)
        if nbad gt 0 then $
           oplot, wave[wbad], check.indiv[iobj,wbad], thick=2, ps=4, color=red
     endfor
     ps_end, /png

     ;; Compute at consistent (within variations of exposure time)
     ;; 2D spec scale
     for iobj=0, nobj-1 do begin
        xpixrange = where(wave ge xrange[0] and wave le xrange[1])
        spec = *objstruct[iobj]
        col = median((spec.xpos)[xpixrange])
        specthumb = transpose((*scistruct[iobj])[col-40:col+39, $
                                                 min(xpixrange):max(xpixrange)])
        if n_elements(specthumbs) eq 0 then begin
           specthumbs = specthumb
        endif else begin
           specthumbs = [specthumbs,specthumb]
        endelse
     endfor
     minthumbs = min(specthumbs, /nan)
     zrangespec = percentile(specthumbs-minthumbs, [0.05, 0.95])
     zrangespec[1] += 0.3*(zrangespec[1]-zrangespec[0])

     ;; Series of plots with coadd and indivs one at a time
     for iobj=0, nobj-1 do begin
        pos1=[0.12, 0.10, 0.97, 0.60]
        pos2=[0.12, 0.605, 0.97, 0.75]
        pos3=[0.12, 0.755, 0.97, 0.90]
        ps_start, filename='LickEPS/'+SNoutname+'.'+setup+'.' $
                  +lname[i]+'.exp'+strn(iobj,f='(I02)')+'.eps'
        plot, [0], /nodata, xrange=xrange, yrange=yrange, $
              xstyle=5, ystyle=1, position=pos1
        xyouts, 0.5, 0.95, SNoutname, align=0.5, /normal
        xyouts, 0.12, 0.95, textoidl(lnametex[i]), align=0, /normal
        xyouts, 0.97, 0.95, setup+'/'+repstr(repstr(sxpar(h, 'INFILE'+strn(iobj,f='(I02)')), $
                                                    '.newwave', ''), $
                                             'OSIRIS-OsirisLongSlitSpectroscopy', ''), $
                align=1, /normal, charsize=1.0
        polyfill, [bc0[i], bc0[i], bc1[i], bc1[i], bc0[i]], $
                  [yrange[0], yrange[1], yrange[1], yrange[0], yrange[0]], $
                  color=skyblue
        polyfill, [rc0[i], rc0[i], rc1[i], rc1[i], rc0[i]], $
                  [yrange[0], yrange[1], yrange[1], yrange[0], yrange[0]], $
                  color=skyblue
        polyfill, [bp0[i], bp0[i], bp1[i], bp1[i], bp0[i]], $
                  [yrange[0], yrange[1], yrange[1], yrange[0], yrange[0]], $
                  color=papaya
        axis, xaxis=0, xrange=xrange, xstyle=1
;        axis, xaxis=1, xrange=xrange*(1.0+galstruct[0].z), xstyle=1
        axis, xaxis=1, xrange=xrange, xstyle=1, xtickname=replicate(' ',30)
        oplot, wave, check.coadd, thick=6, ps=10
        w1 = where(check.indiv eq 0.0) ;; don't plot fluxes of 0.0
        if w1[0] ne -1 then check.indiv[w1] = !values.f_nan
        oplot, wave, check.indiv[iobj,*], thick=2, ps=10, color=darkgreen
        wbad = where(check.mask[iobj,*], nbad)
        if nbad gt 0 then $
           oplot, wave[wbad], check.indiv[iobj,wbad], thick=2, ps=4, color=red

        ;; find thumbnail image boundaries
        xpixrange = where(wave ge xrange[0] and wave le xrange[1])
        spec = *objstruct[iobj]
        col = median((spec.xpos)[xpixrange])
        if min(wave) gt xrange[0] then begin
           pos2[0] = 0.12+(min(wave)-xrange[0])/(xrange[1]-xrange[0])*(0.97-0.12)
           pos3[0] = pos2[0]
        endif
        if max(wave) lt xrange[1] then begin
           pos2[2] = 0.12+(max(wave)-xrange[0])/(xrange[1]-xrange[0])*(0.97-0.12)
           pos3[2] = pos2[2]
        endif

        ;; add pixflat
        pixthumb = transpose(pixflat[col-40:col+39, min(xpixrange):max(xpixrange)])
        zrange = [0.9, 1.1]
        plot, [0], /nodata, /noerase, $
              xrange=xrange, yrange=[-40,40], xstyle=5, ystyle=5, $
              position=pos2
        tvimage, bytscl(pixthumb, min=zrange[0], max=zrange[1], top=251), $
                 position=pos2

        ;; add 2d spec
        specthumb = transpose((*scistruct[iobj])[col-40:col+39, $
                                                 min(xpixrange):max(xpixrange)])
        plot, [0], /nodata, /noerase, $
              xrange=xrange, yrange=[-40,40], xstyle=5, ystyle=5, $
              position=pos3
;        tvimage, bytscl(specthumb-minthumbs, min=zrangespec[0], max=zrangespec[1], top=251), $
;                 position=pos3
        tvimage, bytscl((specthumb-minthumbs)^(0.3), min=(zrangespec[0])^(0.3), $
                        max=(zrangespec[1])^(0.3), top=251), $
                 position=pos3

        ps_end, /png
     endfor
  endfor
  set_plot, 'X'
end
