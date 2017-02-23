;+
; NAME:
;   build_starflux
;
; PURPOSE:
;   Load supernovae and galaxy templates for fitting spectra.
;
; CALLING SEQUENCE:
;   build_starflux, starflux, starwave=starwave, hyper=hyper, sn1a=sn1a, $
;                   sn1bc=sn1bc, sn2l=sn2l, sn2n=sn2n, sn2p=sn2p, $
;                   sn91bg=sn91bg, sn91t=sn91t, eigengal=eigengal, $
;                   eigendir=eigendir, gal=gal, sne=sne, clobber=clobber, $
;                   silent=silent, allhyper=allhyper, allsn1a=allsn1a, $
;                   allsn2l=allsn2l, allsn2n=allsn2n, allsn2p=allsn2p, $
;                   allsn91bg=allsn91bg, allsn91t=allsn91t, allsne=allsne, $
;                   catalog=catalog
;
; INPUTS:
;
; REQUIRED KEYWORDS:
;
; OPTIONAL KEYWORDS:
;   starwave   - if initially empty, then starwave is set to the
;                wavelength range present in D. Schlegel's PCA galaxy
;                templates.  If initially non-empty, then the output
;                templates are interpolated onto this wavelength
;                range.
;   clobber    - Instead of appending output to starflux/catalog,
;                overwrite them.
;   silent     - Silence output to terminal
;   hyper,sn1a,sn1bc,sn2l,sn2n,sn2p,sn91bg,sn91t
;              - request Nugent02 template for specific epoch (in
;                SN-frame days past max)
;   eigengal   - set to filename containing D. Schlegel's PCA galaxy templates.
;   eigendir   - set directory for above.
;   gal        - request D. Schlegel's PCA galaxy templates.
;   sne        - request "standard" supernovae templates: sn1a=20,
;                sn1a=65, sn1bc=15, sn2p=61
;   allhyper, allsn1a, allsn2l, allsn2n, allsn2p, allsn91bg, allsn91t
;              - request all epochs of given Nugent02 template series
;   allsne     - request all epochs of all Nugent02 template series's.
;   
; OUTPUTS:
;   starflux   - array onto which requested templates are appended.
;
; OPTIONAL OUTPUTS:
;   catalog    - string array onto which requested templates
;                descriptions are appended.
;   
;
; COMMENTS:
;   Output is normalized in the unit-vector sense.
; EXAMPLES:
;
; BUGS:
;   Should update to include Hsiao templates.
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   05-15-07 - Written by Joshua E. Meyers - LBNL
;------------------------------------------------------------------------------
pro build_starflux, starflux, starwave=starwave, hyper=hyper, sn1a=sn1a, $
                    sn1bc=sn1bc, sn2l=sn2l, sn2n=sn2n, sn2p=sn2p, $
                    sn91bg=sn91bg, sn91t=sn91t, eigengal=eigengal, $
                    eigendir=eigendir, gal=gal, sne=sne, clobber=clobber, $
                    silent=silent, allhyper=allhyper, allsn1a=allsn1a, $
                    allsn1bc=allsn1bc, allsn2l=allsn2l, allsn2n=allsn2n, $
                    allsn2p=allsn2p, allsn91bg=allsn91bg, allsn91t=allsn91t, $
                    allsne=allsne, catalog=catalog

  ;----------
  ; Start from scratch if clobber is set; otherwise append results.

  if (keyword_set(clobber)) then begin
     delvarx, starflux
     delvarx, starwave
     delvarx, catalog
  end

  ;----------
  ; Read in SNe templates
  
  dir=concat_dir(getenv('JEM_REF'),'sn_templates')
  file='sne-templates.sav'
  allfiles=file_search(djs_filepath(file, root_dir = dir), count=ct)
  if (ct eq 0) then $
     message, 'Unable to find SNEFILE matching '+file
  thisfile=allfiles[(reverse(sort(allfiles)))[0]]
  if (~keyword_set(silent)) then $
     print, 'Selecting SNe template file=' + thisfile
  restore, thisfile

  ;----------
  ; Default galaxy template.

  if (n_elements(eigendir) eq 0) then $
     eigendir = concat_dir(getenv('IDLSPEC2D_DIR'),'templates')
  if (n_elements(eigengal) eq 0) then $
     eigengal = 'spEigenGal-*.fits'

  ;----------
  ; Use the wavelengths given by galaxy file unless user supplies starwave.

  if (NOT keyword_set(starwave)) then begin
     allfiles=file_search(djs_filepath(eigengal, root_dir = eigendir), $
                          count=ct)
     if (ct eq 0) then $
        message, 'Unable to find EIGENFILE matching '+eigengal
     galfile=allfiles[(reverse(sort(allfiles)))[0]]
     if (~keyword_set(silent)) then $
        print, 'Selecting EIGENFILE=' + galfile
     s0=readfits(galfile, shdr, /silent)
     starloglam0=sxpar(shdr, 'COEFF0')
     stardloglam=sxpar(shdr, 'COEFF1')
     nstarpix=(size(s0,/dimens))[0]
     starwave=10^(dindgen(nstarpix)*stardloglam+starloglam0)
  endif

  ;----------
  ; Add galaxy spectra.  Use already opened galaxy file if possible.

  if (keyword_set(gal)) then begin
     if (NOT keyword_set(galfile)) then begin
        allfiles=file_search(djs_filepath(eigengal, root_dir=eigendir), $
                             count=ct)
        if (ct eq 0) then $
           message, 'Unable to find EIGENFILE matching '+eigengal
        galfile=allfiles[(reverse(sort(allfiles)))[0]]
        if (~keyword_set(silent)) then $
           print, 'Selecting EIGENFILE='+galfile
        ss0=readfits(galfile, shdr, /silent)
        starloglam0=sxpar(shdr, 'COEFF0')
        stardloglam=sxpar(shdr, 'COEFF1')
        nstarpix=(size(ss0,/dimens))[0]
        if(size(ss0,/n_dimen) eq 1) then nstar = 1 $
        else nstar=(size(ss0,/dimens))[1]
        swave=10^(dindgen(nstarpix)*stardloglam+starloglam0)
        ; need to interpolate from swave -> starwave
        s0=dblarr(n_elements(starwave),nstar)
        for istar=0,nstar-1 do $
           s0[*,istar]=interpol(ss0[*,istar], swave, starwave)
     endif
     if (keyword_set(starflux)) then $
        starflux = [[starflux],[s0]] $
     else $
        starflux = s0
     if (keyword_set(catalog)) then $
        catalog=[catalog,['PCAgal1','PCAgal2','PCAgal3','PGAgal4']] $
     else $
        catalog=['PCAgal1','PCAgal2','PCAgal3','PGAgal4']
  endif

  ;----------
  ; Add 'standard' SNe spectra
  
  if (keyword_set(sne)) then begin
     build_starflux, starflux, starwave=starwave, sn1a=20, catalog=catalog, $
                     silent=silent
     build_starflux, starflux, starwave=starwave, sn1a=65, sn1bc=15, sn2p=61, $
       catalog=catalog, silent=silent
  endif

  ;----------
  ; Add hypernova spectrum.

  if (n_elements(hyper) ne 0) then begin
     if (~keyword_set(hday)) then restore, savfile
     if (hyper ge max(hday)) then begin
        if (~keyword_set(silent)) then $
           print, 'hyper > max(days), using max = ', max(hday)
        d=max(hday)
     endif else if (hyper le min(hday)) then begin
        if (~keyword_set(silent)) then $
           print, 'hyper < min(days), using min = ', min(hday)
        d=min(hday)
     endif else begin
        ind=0
        while (hday[ind] lt hyper) do ind+=1
        if (hday[ind]-hyper lt hyper-hday[ind-1]) then d=hday[ind] $
        else d=hday[ind-1]
        if (~keyword_set(silent)) then $
           print, 'using d = ', d, ' for hyper'
     endelse
     i=where(hday eq d)*n_elements(wave)
     sflux=hflux[i:i+n_elements(wave)-1]
     s0=interpol(sflux,wave,starwave)
     if (keyword_set(starflux)) then $
        starflux = [[starflux],[s0]] $
     else $
        starflux = s0
     if (keyword_set(catalog)) then $
        catalog=[catalog,'hyper'+string(d)] $
     else $
        catalog='hyper'+string(d)
  endif

  ;----------
  ; Add all hypernova spectra
  if (keyword_set(allhyper)) then begin
     if (~keyword_set(hday)) then restore, savfile
     for id=0,n_elements(hday)-1 do begin
        build_starflux,starflux,starwave=starwave,hyper=hday[id], $
                       silent=silent, catalog=catalog
     endfor
  endif


  ;----------
  ; Add SN type Ia spectrum.

  if (n_elements(sn1a) ne 0) then begin
     if (~keyword_set(sn1aday)) then restore, savfile
     if (sn1a ge max(sn1aday)) then begin
        if (~keyword_set(silent)) then $
           print, 'sn1a > max(days), using max = ', max(sn1aday)
        d=max(sn1aday)
     endif else if (sn1a le min(sn1aday)) then begin
        if (~keyword_set(silent)) then $
           print, 'sn1a < min(days), using min = ', min(sn1aday)
        d=min(sn1aday)
     endif else begin
        ind=0
        while (sn1aday[ind] lt sn1a) do ind+=1
        if (sn1aday[ind]-sn1a lt sn1a-sn1aday[ind-1]) then d=sn1aday[ind] $
        else d=sn1aday[ind-1]
        if (~keyword_set(silent)) then $
           print, 'using d = ', d, ' for sn1a'
     endelse
     i=where(sn1aday eq d)*n_elements(wave)
     sflux=sn1aflux[i:i+n_elements(wave)-1]
     s0=interpol(sflux,wave,starwave)
     if (keyword_set(starflux)) then $
        starflux = [[starflux],[s0]] $
     else $
        starflux = s0
     if (keyword_set(catalog)) then $
        catalog=[catalog,'sn1a'+string(d)] $
     else $
        catalog='sn1a'+string(d)
  endif

  ;----------
  ; Add all SN type Ia spectra
  if (keyword_set(allsn1a)) then begin
     if (~keyword_set(sn1aday)) then restore, savfile
     for id=0,n_elements(sn1aday)-1 do begin
        build_starflux,starflux,starwave=starwave,sn1a=sn1aday[id], $
                       catalog=catalog, silent=silent
    endfor
  endif


  ;----------
  ; Add SN type Ibc spectrum.

  if (n_elements(sn1bc) ne 0) then begin
     if (~keyword_set(sn1bcday)) then restore, savfile
     if (sn1bc ge max(sn1bcday)) then begin
        if (~keyword_set(silent)) then $
           print, 'sn1bc > max(days), using max = ', max(sn1bcday)
        d=max(sn1bcday)
     endif else if (sn1bc le min(sn1bcday)) then begin
        if (~keyword_set(silent)) then $
           print, 'sn1bc < min(days), using min = ', min(sn1bcday)
        d=min(sn1bcday)
     endif else begin
        ind=0
        while (sn1bcday[ind] lt sn1bc) do ind+=1
        if (sn1bcday[ind]-sn1bc lt sn1bc-sn1bcday[ind-1]) then d=sn1bcday[ind] $
        else d=sn1bcday[ind-1]
        if (~keyword_set(silent)) then $
           print, 'using d = ', d, ' for sn1bc'
     endelse
     i=where(sn1bcday eq d)*n_elements(wave)
     sflux=sn1bcflux[i:i+n_elements(wave)-1]
     s0=interpol(sflux,wave,starwave)
     if (keyword_set(starflux)) then $
        starflux = [[starflux],[s0]] $
     else $
        starflux = s0
     if (keyword_set(catalog)) then $
        catalog=[catalog,'sn1bc'+string(d)] $
     else $
        catalog='sn1bc'+string(d)
  endif

  ;----------
  ; Add all SN type Ibc spectra
  if (keyword_set(allsn1bc)) then begin
     if (~keyword_set(sn1bcday)) then restore, savfile
     for id=0,n_elements(sn1bcday)-1 do begin
        build_starflux,starflux,starwave=starwave,sn1bc=sn1bcday[id], $
                       catalog=catalog, silent=silent
     endfor
  endif


  ;----------
  ; Add SN type IIl spectrum.

  if (n_elements(sn2l) ne 0) then begin
     if (~keyword_set(sn2lday)) then restore, savfile
     if (sn2l ge max(sn2lday)) then begin
        if (~keyword_set(silent)) then $
           print, 'sn2l > max(days), using max = ', max(sn2lday)
        d=max(sn2lday)
     endif else if (sn2l le min(sn2lday)) then begin
        if (~keyword_set(silent)) then $
           print, 'sn2l < min(days), using min = ', min(sn2lday)
        d=min(sn2lday)
     endif else begin
        ind=0
        while (sn2lday[ind] lt sn2l) do ind+=1
        if (sn2lday[ind]-sn2l lt sn2l-sn2lday[ind-1]) then d=sn2lday[ind] $
        else d=sn2lday[ind-1]
        if (~keyword_set(silent)) then $
           print, 'using d = ', d, ' for sn2l'
     endelse
     i=where(sn2lday eq d)*n_elements(wave)
     sflux=sn2lflux[i:i+n_elements(wave)-1]
     s0=interpol(sflux,wave,starwave)
     if (keyword_set(starflux)) then $
        starflux = [[starflux],[s0]] $
     else $
        starflux = s0
     if (keyword_set(catalog)) then $
        catalog=[catalog,'sn2l'+string(d)] $
     else $
        catalog='sn2l'+string(d)
  endif

  ;----------
  ; Add all SN type IIl spectra
  if (keyword_set(allsn2l)) then begin
     if (~keyword_set(sn2lday)) then restore, savfile
     for id=0,n_elements(sn2lday)-1 do begin
        build_starflux,starflux,starwave=starwave,sn2l=sn2lday[id], $
                       catalog=catalog,silent=silent
     endfor
  endif


  ;----------
  ; Add SN type IIn spectrum.

  if (n_elements(sn2n) ne 0) then begin
     if (~keyword_set(sn2nday)) then restore, savfile
     if (sn2n ge max(sn2nday)) then begin
        if (~keyword_set(silent)) then $
           print, 'sn2n > max(days), using max = ', max(sn2nday)
        d=max(sn2nday)
     endif else if (sn2n le min(sn2nday)) then begin
        if (~keyword_set(silent)) then $
           print, 'sn2n < min(days), using min = ', min(sn2nday)
        d=min(sn2nday)
     endif else begin
        ind=0
        while (sn2nday[ind] lt sn2n) do ind+=1
        if (sn2nday[ind]-sn2n lt sn2n-sn2nday[ind-1]) then d=sn2nday[ind] $
        else d=sn2nday[ind-1]
        if (~keyword_set(silent)) then $
           print, 'using d = ', d, ' for sn2n'
     endelse
     i=where(sn2nday eq d)*n_elements(wave)
     sflux=sn2nflux[i:i+n_elements(wave)-1]
     s0=interpol(sflux,wave,starwave)
     if (keyword_set(starflux)) then $
        starflux = [[starflux],[s0]] $
     else $
        starflux = s0
     if (keyword_set(catalog)) then $
        catalog=[catalog,'sn2n'+string(d)] $
     else $
        catalog='sn2n'+string(d)
  endif

  ;----------
  ; Add all SN type IIn spectra
  if (keyword_set(allsn2n)) then begin
     if (~keyword_set(sn2nday)) then restore, savfile
     for id=0,n_elements(sn2nday)-1 do begin
        build_starflux,starflux,starwave=starwave,sn2n=sn2nday[id], $
                       catalog=catalog, silent=silent
     endfor
  endif


  ;----------
  ; Add SN type IIp spectrum.

  if (n_elements(sn2p) ne 0) then begin
     if (~keyword_set(sn2pday)) then restore, savfile
     if (sn2p ge max(sn2pday)) then begin
        if (~keyword_set(silent)) then $
           print, 'sn2p > max(days), using max = ', max(sn2pday)
        d=max(sn2pday)
     endif else if (sn2p le min(sn2pday)) then begin
        if (~keyword_set(silent)) then $
           print, 'sn2p < min(days), using min = ', min(sn2pday)
        d=min(sn2pday)
     endif else begin
        ind=0
        while (sn2pday[ind] lt sn2p) do ind+=1
        if (sn2pday[ind]-sn2p lt sn2p-sn2pday[ind-1]) then d=sn2pday[ind] $
        else d=sn2pday[ind-1]
        if (~keyword_set(silent)) then $
           print, 'using d = ', d, ' for sn2p'
     endelse
     i=where(sn2pday eq d)*n_elements(wave)
     sflux=sn2pflux[i:i+n_elements(wave)-1]
     s0=interpol(sflux,wave,starwave)
     if (keyword_set(starflux)) then $
        starflux = [[starflux],[s0]] $
     else $
        starflux = s0
     if (keyword_set(catalog)) then $
        catalog=[catalog,'sn2p'+string(d)] $
     else $
        catalog='sn2p'+string(d)
  endif

  ;----------
  ; Add all SN type IIp spectra
  if (keyword_set(allsn2p)) then begin
     if (~keyword_set(sn2pday)) then restore, savfile
     for id=0,n_elements(sn2pday)-1 do begin
        build_starflux,starflux,starwave=starwave,sn2p=sn2pday[id], $
                       catalog=catalog,silent=silent
     endfor
  endif


  ;----------
  ; Add SN type 91bg spectrum.

  if (n_elements(sn91bg) ne 0) then begin
     if (~keyword_set(sn91bgday)) then restore, savfile
     if (sn91bg ge max(sn91bgday)) then begin
        if (~keyword_set(silent)) then $
           print, 'sn91bg > max(days), using max = ', max(sn91bgday)
        d=max(sn91bgday)
     endif else if (sn91bg le min(sn91bgday)) then begin
        if (~keyword_set(silent)) then $
           print, 'sn91bg < min(days), using min = ', min(sn91bgday)
        d=min(sn91bgday)
     endif else begin
        ind=0
        while (sn91bgday[ind] lt sn91bg) do ind+=1
        if (sn91bgday[ind]-sn91bg lt sn91bg-sn91bgday[ind-1]) then $
           d=sn91bgday[ind] $
        else d=sn91bgday[ind-1]
        if (~keyword_set(silent)) then $
           print, 'using d = ', d, ' for sn91bg'
     endelse
     i=where(sn91bgday eq d)*n_elements(wave)
     sflux=sn91bgflux[i:i+n_elements(wave)-1]
     s0=interpol(sflux,wave,starwave)
     if (keyword_set(starflux)) then $
        starflux = [[starflux],[s0]] $
     else $
        starflux = s0
     if (keyword_set(catalog)) then $
        catalog=[catalog,'sn91bg'+string(d)] $
     else $
        catalog='sn91bg'+string(d)
  endif

  ;----------
  ; Add all SN type 91bg spectra
  if (keyword_set(allsn91bg)) then begin
     if (~keyword_set(sn91bgday)) then restore, savfile
     for id=0,n_elements(sn91bgday)-1 do begin
        build_starflux,starflux,starwave=starwave,sn91bg=sn91bgday[id], $
                       catalog=catalog, silent=silent
     endfor
  endif


  ;----------
  ; Add SN type 91t spectrum.

  if (n_elements(sn91t) ne 0) then begin
     if (~keyword_set(sn91tday)) then restore, savfile
     if (sn91t ge max(sn91tday)) then begin
        if (~keyword_set(silent)) then $
           print, 'sn91t > max(days), using max = ', max(sn91tday)
        d=max(sn91tday)
     endif else if (sn91t le min(sn91tday)) then begin
        if (~keyword_set(silent)) then $
           print, 'sn91t < min(days), using min = ', min(sn91tday)
        d=min(sn91tday)
     endif else begin
        ind=0
        while (sn91tday[ind] lt sn91t) do ind+=1
        if (sn91tday[ind]-sn91t lt sn91t-sn91tday[ind-1]) then d=sn91tday[ind] $
        else d=sn91tday[ind-1]
        if (~keyword_set(silent)) then $
           print, 'using d = ', d, ' for sn91t'
     endelse
     i=where(sn91tday eq d)*n_elements(wave)
     sflux=sn91tflux[i:i+n_elements(wave)-1]
     s0=interpol(sflux,wave,starwave)
     if (keyword_set(starflux)) then $
        starflux = [[starflux],[s0]] $
     else $
        starflux = s0
     if (keyword_set(catalog)) then $
        catalog=[catalog,'sn91t'+string(d)] $
     else $
        catalog='sn91t'+string(d)
  endif

  ;----------
  ; Add all SN type 91t spectra
  if (keyword_set(allsn91t)) then begin
     if (~keyword_set(sn91tday)) then restore, savfile
     for id=0,n_elements(sn91tday)-1 do begin
        build_starflux,starflux,starwave=starwave,sn91t=sn91tday[id], $
                       catalog=catalog,silent=silent
     endfor
  endif


  ;----------
  ; Add all SN spectra!
  if (keyword_set(allsne)) then begin
     build_starflux, starflux, starwave=starwave, /allhyper, /allsn1a, $
       /allsn1bc, /allsn2p, /allsn2l, /allsn2n, /allsn91bg, /allsn91t, $
       catalog=catalog, silent=silent
  endif

  ;----------
  ; Normalize starflux vectors.

  ndim=(size(starflux,/n_dimen))
  if (ndim eq 1) then nstar=1 $
  else nstar=(size(starflux,/dimens))[1]
  for istar=0,nstar-1 do $
     starflux[*,istar] /= sqrt(total(starflux[*,istar]*starflux[*,istar]))
end
