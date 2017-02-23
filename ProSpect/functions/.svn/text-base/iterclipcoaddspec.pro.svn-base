Function IterClipCoaddSpec, spectra, $
                            frame=frame, $
                            verbose=verbose, $
                            norm=norm, $
                            dlam=dlam, $
                            dloglam=dloglam

  waveminmax=[1.e19, -1.e19]
  for ispec=0, n_elements(spectra)-1 do begin
     if size( spectra[ispec], /tname ) eq 'POINTER' then begin
        if ~ptr_valid(spectra[ispec]) then continue
        waveminmax1=minmax((*spectra[ispec])->wavelength(frame=frame))
     endif else begin
        if ~obj_valid(spectra[ispec]) then continue
        waveminmax1=minmax((spectra[ispec])->wavelength(frame=frame))
     endelse
     waveminmax[0] = waveminmax[0] < waveminmax1[0]
     waveminmax[1] = waveminmax[1] > waveminmax1[1]
  endfor
  
  print, waveminmax

  if n_elements(dlam) ne 0 then begin
     nbins = (waveminmax[1]-waveminmax[0])/dlam+1
     binbounds = findgen(nbins)*dlam+waveminmax[0]
  endif else begin
     nbins = ((alog10(waveminmax[1])-alog10(waveminmax[0]))/dloglam)+1
     binbounds = 10^(findgen(nbins)*dloglam+alog10(waveminmax[0]))
  endelse
  
  print, minmax(binbounds)

  for ispec=0, n_elements(spectra)-1 do begin
     if keyword_set(verbose) then print, 'spectrum'+string(ispec+1, format='(I)')
     if size( spectra[ispec], /tname ) eq 'POINTER' then begin
        if ~ptr_valid(spectra[ispec]) then continue
        wave=(*spectra[ispec])->wavelength(frame=frame)
        flux=(*spectra[ispec])->flux(frame=frame)
        ivar=(*spectra[ispec])->ivar(frame=frame)
     endif else begin
        if ~obj_valid(spectra[ispec]) then continue
        wave=spectra[ispec]->wavelength(frame=frame)
        flux=spectra[ispec]->flux(frame=frame)
        ivar=spectra[ispec]->ivar(frame=frame)
     endelse
     winfinite=where(1-finite(ivar))
     if winfinite[0] ne -1 then ivar[winfinite]=0.
     if n_elements(norm) ne 0 then begin
        wflux=where( wave ge norm[0] $
                     and wave le norm[1] )
        if n_elements(normvalue) eq 0 then $
           normvalue = total(flux[wflux]) $
        else begin
           flux *= normvalue/total(flux[wflux])
           ivar /= (normvalue/total(flux[wflux]))^2
        endelse
     endif
     binflux=dblarr(n_elements(binbounds)-1)
     binivar=dblarr(n_elements(binbounds)-1)
     jem_bin_spectrum, wave, flux, ivar, $
                       binflux=binflux, binivar=binivar, $
                       binbounds=binbounds, /append
     if n_elements(fluxes) eq 0 then begin
        fluxes = binflux 
        ivars = binivar
     endif else begin
        fluxes = [[[fluxes]], [[binflux]]]
        ivars  = [[[ivars]], [[binivar]]]
     endelse
  endfor
  spec = avsigclip( fluxes, ivars, 3. )
  binwave = 0.5*(shift(binbounds,1)+binbounds)[1:n_elements(binbounds)-1]
  outspec = obj_new('spectrum', binwave, spec.flux, ivar=spec.ivar)
  return, outspec
end




;  function avsigclip, array, invvar, threshold, sigmap, $
;                    inmask=inmask,niteration=niteration
