;+
; NAME:
;   jem_zcompute2
;
; PURPOSE:
;   Find possible redshift matches for a set of spectra using a set of
;   eigen-templates.
;
; CALLING SEQUENCE
;   result = jem_zcompute2( objwave, objflux, objivar, $
;      starflux=starflux, logwave0=logwave0, dlogwave=dlogwave, $
;      zmin=zmin, zmax=zmax, [ npoly=npoly, nfind=nfind, $
;      mindof=mindof, chi2arr=chi2arr, dofarr=dofarr, ztry=ztry, $
;      silent=silent, synflux=synflux ]
;
Function create_zans, nstar, nfind
  zans1 = create_struct( $
          name = 'ZANS'+strtrim(string(nstar),1), $
          'z', 0.0, $
          'z_err', 0.0, $
          'chi2', 0.0, $
          'dof', 0L, $
          'theta', fltarr(nstar), $
          'theta_covar', fltarr(nstar,nstar))
  return, replicate(zans1,nfind)
end
;------------------------------------------------------------------------------------
Function jem_zcompute2, objwave, objflux, objivar, $
                        starflux=starflux, logwave0=logwave0, dlogwave=dlogwave, $
                        zmin=zmin, zmax=zmax, npoly=npoly, nfind=nfind, $
                        minsep=minsep, width=width, mindof=mindof, $
                        chi2arr=chi2arr, dofarr=dofarr, ztry=ztry, $
                        silent=silent, synflux=synflux

  if n_elements(nfind) eq 0 then nfind=1
  if n_elements(mindof) eq 0 then mindof=20
  if n_elements(npoly) eq 0 then npoly=0
  if n_elements(width) eq 0 then width=0.01
  if n_elements(minsep) eq 0 then minsep=width

  ;---------
  ; Get dimensions

  ndim=size( objflux, /n_dimen )
  if (ndim eq 1) then nobj=1 $
  else nobj=(size( objflux, /dimen ))[1]
  nobjpix=(size( objflux, /dimen ))[0]
  ndim=size( starflux, /n_dimen )
  if (ndim eq 1) then nstar=1 $
  else nstar=(size( starflux, /dimen ))[1]
  nstarpix=(size( starflux, /dimen))[0]

  if n_elements(objivar) eq 0 then objivar=dblarr( nobjpix, nobj )+1.0

  if (~keyword_set(silent)) then begin
     t0 = systime(1)
     print, 'Redshifting '+strtrim( string(nobj), 2 )+' spectra.'
  endif

  ;-------------
  ; Rebin spectra into log-linear wavelengths
  print, 'Rebinning spectra into log-linear wavelengths'
  lmin = min(objwave)
  lmax = max(objwave)
  nbin = ceil((alog10(lmax/lmin)/dlogwave))+1
  binwave = 10.^(findgen(nbin)*dlogwave+alog10(lmin))
  objbins = wave2bins(binwave)
  binflux = dblarr( nbin, nobj )
  binivar = dblarr( nbin, nobj )
  for iobj=0, nobj-1 do begin
     objflux1=objflux[*,iobj]
     objivar1=objivar[*,iobj]
;     jem_drizzle1d, objwave, objflux1, ivar=objivar1, $
;                    outflux=outflux, outivar=outivar, $
;                    outbins=objbins
     drizzle1d, objwave, objflux1, ivar=objivar1, $
                outflux=outflux, outivar=outivar, $
                outbins=objbins
     binflux[*,iobj]=outflux
     binivar[*,iobj]=outivar
  endfor

  ;----------
  ; Create master ztry table
  print, 'Creating master ztry table'
  templatewave = 10.^(dindgen(nstarpix)*dlogwave+logwave0)
  zminposs = (binwave[mindof-1]/templatewave[nstarpix-1])-1.
  zmaxposs = (binwave[nbin-1]/templatewave[mindof-1])-1.
  nzposs = ceil(alog10((1.+zmaxposs)/(1.+zminposs))/dlogwave)+1
  zpossible = 10.^(dindgen(nzposs)*dlogwave+alog10(zminposs+1))-1
  ztry=zpossible
  if n_elements(zmin) ne 0 then ztry=ztry[where(ztry ge min(zmin))] $
  else zmin=fltarr(nobjs)+min(ztry)
  if n_elements(zmax) ne 0 then ztry=ztry[where(ztry le max(zmax))] $
  else zmax=fltarr(nobjs)+max(ztry)
  nz=n_elements(ztry)

  ;------------
  ; Determine min/max z-indices for objects
  print, 'Determining min/max z-indices'
  minind=intarr(nobj)
  maxind=intarr(nobj)
  for iobj=0, nobj-1 do begin
     minind[iobj]=min(where(ztry ge zmin[iobj]))
     maxind[iobj]=max(where(ztry le zmax[iobj]))
  endfor
  if ~keyword_set(silent) then begin
     print, 'zmin: '+strtrim( string(min(zmin)), 2 )
     print, 'zmax: '+strtrim( string(max(zmax)), 2 )
  endif

  ;-------------
  ; Add Chebyshev polynomials
  print, 'Adding Chebyshev polynomials
  if (npoly ne 0) then begin
     fixed_template = cheb_array( npoly-1, alog10(templatewave) )
     fixed_template2 = cheb_array( npoly-1, objwave )
  endif

  ;------------
  ; Initialize arrays
  print, 'Initializing arrays'
  chi2arr = fltarr( nz, nobj )
  dofarr  = fltarr( nz, nobj )
  synflux = fltarr( nobjpix, nfind, nobj )
  sqobjivar = sqrt( binivar )
  starfluxz = fltarr( nobjpix, nstar )
  for ishift=0, nz-1 do begin
     z=ztry[ishift]
     counter, ishift, nz-1, 'z='+string(ztry[ishift],f='(F9.3)')+'  '
     shift=round(alog10((1+z)/(binwave[0]/templatewave[0]))/dlogwave)
     template1 = shift(starflux,shift,0)
     if shift le 0 then begin
        minpix = 0
        maxpix = (nstarpix+shift-1) < (nbin-1)
     endif else begin
        minpix = shift
        maxpix = (nstarpix-1) < (nbin-1)
     endelse
     for iobj=0L, nobj-1 do begin
        if ishift ge minind[iobj] and ishift le maxind[iobj] then begin
           idx = indgen(maxpix-minpix+1)+minpix
           if keyword_set(fixed_template) then $
              chi2arr[ishift,iobj] = $
              computechi2(binflux[idx,iobj],sqobjivar[idx,iobj], $
                          ([[template1],[fixed_template]])[idx,*], $
                          dof=dof) $
           else $
              chi2arr[ishift,iobj] = $
              computechi2(binflux[idx,iobj],sqobjivar[idx,iobj], $
                          template1[idx,*], dof=dof)
           dofarr[ishift,iobj]=dof
        endif
     endfor
  endfor
  for iobj=0L, nobj-1 do begin
     zans1=create_zans(nstar+npoly,nfind)
     indx = where(dofarr[*,iobj] ge mindof, ct)
     if ct ge width then begin
        xpeak1 = find_nminima(chi2arr[indx,iobj], ztry[indx], $
                              dofarr=dofarr[indx,iobj], nfind=nfind, $
                              minsep=minsep, width=width, ypeak=ypeak1, $
                              xerr=xerr1, npeak=npeak, errcode=errcode)
        zans1[0:npeak-1].z=xpeak1
        zans1[0:npeak-1].z_err=xerr1 * (errcode eq 0) + errcode
        zans1[0:npeak-1].chi2 = ypeak1
        for ipeak=0L, npeak-1 do begin
           junk=min(abs(ztry-xpeak1[ipeak]),iz)
           zans1[ipeak].dof=dofarr[iz,iobj]
        endfor
        zans1[0:npeak-1].chi2 *= zans1[0:npeak-1].dof
     endif else if (ct ge 1) then begin
        npeak=1
        zans1[0].chi2 = -min(-chi2arr[indx]/dofarr[indx],iz)
        zans1[0].z = ztry[indx[iz]]
        zans1[0].z_err = 0
        zans1[0].dof = dofarr[indx[iz],iobj]
        zans1[0].chi2 *= zans1[0].dof
     endif else begin
        npeak=0
     endelse

     for ipeak=0L,npeak-1 do begin
        z=xpeak1[ipeak]
        starwavez=templatewave*(1.0+z)
        for istar=0, nstar-1 do $
           starfluxz[*,istar] = interpol(starflux[*,istar], starwavez, objwave)

        idx = where(objwave le max(starwavez) and objwave ge min(starwavez))
        if keyword_set(fixed_template2) then begin
           thischi2 = computechi2(objflux[idx,iobj],sqrt(objivar[idx,iobj]), $
                                  ([[starfluxz],[fixed_template2]])[idx,*], acoeff=acoeff, $
                                  covar=covar, yfit=synflux1)
           zans1[ipeak].theta=acoeff
           zans1[ipeak].theta_covar=covar
           synflux[idx,ipeak,iobj]=synflux1
        endif else begin
           thischi2 = computechi2(objflux[idx,iobj],sqrt(objivar[idx,iobj]), $
                                  starfluxz[idx,*], acoeff=acoeff, $
                                  covar=covar, yfit=synflux1)
           zans1[ipeak].theta=acoeff
           zans1[ipeak].theta_covar=covar
           synflux[idx,ipeak,iobj]=synflux1
        endelse
     endfor ;end peak loop

     if (iobj eq 0) then zans=zans1 $
     else zans=[zans,zans1]
  endfor
  return, zans
end
