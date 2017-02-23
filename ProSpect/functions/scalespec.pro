Function scalespec, $
   spec, $  ;;spectrum to scale.  will not be overwritten
   refspec=refspec, $ ;;scale to this spectrum
   region=region, $ ;;scale to refspec in this wavelength region
   const=const, $ ;;scale by this constant
   flux=flux ;;scale by flux instead of minimizing chi2
  
  fluxin = spec->flux()
  wavein = spec->wavelength()
  ivarin = spec->ivar()

  if n_elements(const) ne 0 then begin
     ratio = const
  endif else begin
     
     refflux = refspec->flux()
     refwave = refspec->wavelength()
     refivar = refspec->ivar()
     
     if n_elements(region) eq 0 then begin
        region = dblarr(2)
        region[0] = max([min(wavein),min(refwave)])
        region[1] = min([max(wavein),max(refwave)])
     endif
     if keyword_set(flux) then begin
        wref = where(refwave ge region[0] and refwave le region[1])
        ws = where(wavein ge region[0] and refwave le region[1])
        if wref[0] eq -1 or ws[0] eq -1 then $
           message, 'invalid region in scalespec'
        reftotal = total(refflux[wref])
        stotal = total(fluxin[ws])
        ratio = reftotal/stotal
     endif else begin ;;must be chi2 comparison
        outbins = wave2bins(refwave)
        spec1 = rebinspectrum(spec, outbins=outbins)
        flux1 = spec1->flux()
        ivar1 = spec1->ivar()
        ivar2 = 1./(1./(ivar1)+1./(refivar))
        winfinite=where(1-finite(ivar2))
        if winfinite[0] ne -1 then begin
           ivar2[winfinite]=0.
           flux1[winfinite]=0.
        endif
        w=where(refwave ge region[0] and refwave le region[1])
        if w[0] eq -1 then $
           message, 'invalid region in scalespec'
        junk = computechi2(refflux[w], ivar2[w], flux1[w], acoeff=acoeff)
        ratio = acoeff[0]
        obj_destroy, spec1
     endelse
  endelse

  return, obj_new('spectrum', wavein, fluxin*ratio, ivar=ivarin/ratio^2)
end
