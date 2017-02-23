Function masktell, wave
  flag = bytarr(n_elements(wave))
  ;;Telluric lines
  flag or= wave gt 6270.2 and wave lt 6331.7
  flag or= wave gt 6862.1 and wave lt 6964.6 ;; B-band
  flag or= wave gt 7143.3 and wave lt 7398.2 ;;
  flag or= wave gt 7585.8 and wave lt 7703.0 ;; A-band
  flag or= wave gt 7887.6 and wave lt 8045.8
  flag or= wave gt 8916.0 and wave lt 9929.8
  return, 1-flag
end

Function maskstar, wave, starname
  clight = 299792.458d
  flag = bytarr(n_elements(wave))

  case strupcase(strtrim(starname,2)) of
     'FEIGE34' : begin
        ;;Balmer lines
        linewaves = [6564.614, 4862.677, 4341.676, 4102.884]
        ;; some more lines in spectrum of HZ44
        ;linewaves = [linewaves, 5877.071, 6500.0, 6679.669, 7066.928]
     end
     'HZ44' : begin
        ;;Balmer lines
        linewaves = [6564.614, 4862.677, 4341.676, 4102.884]
        ;; some more lines in spectrum of HZ44
        linewaves = [linewaves, 5877.071, 6500.0, 6679.669, 7066.928]
     end
     'BD28' : begin
        ;;Balmer lines
        linewaves = [6564.614, 4862.677, 4341.676, 4102.884]
        ;; some more lines in spectrum of BD28
        linewaves = [linewaves, 5400.0]
     end
     'BD33' : begin
        ;;Balmer lines
        linewaves = [6564.614, 4862.677, 4341.676, 4102.884]
        ;; some more lines in spectrum of BD28
        ;linewaves = [linewaves, 5400.0]
     end
     'BD284211' : begin
        ;;Balmer lines
        linewaves = [6564.614, 4862.677, 4341.676, 4102.884]
        ;; some more lines in spectrum of BD28
        linewaves = [linewaves, 5400.0]
     end
     'LDS749B' : begin
        ;;Balmer lines
        linewaves = [6564.614, 4862.677, 4341.676, 4102.884]
        ;; some more lines in spectrum of BD28
        linewaves = [linewaves, 5880.0, 6675.0]
     end
     'G191B2B' : begin
        ;;Balmer lines
        linewaves = [6564.614, 4862.677, 4341.676, 4102.884]
        ;; some more lines in spectrum of G191B2B
        ;linewaves = [linewaves, 5400.0]
     end
     'GD153' : begin
        ;;Balmer lines
        linewaves = [6564.614, 4862.677, 4341.676, 4102.884]
        ;; some more lines in spectrum of G191B2B
        ;linewaves = [linewaves, 5400.0]
     end
  endcase

  for i=0, n_elements(linewaves)-1 do begin
     flag or= wave gt linewaves[i]*(1d - 1500d/clight) $
              and wave lt linewaves[i]*(1d + 1500d/clight)
  endfor

  ;; clip edges
  flag[0:15] = 1
  flag[n_elements(flag)-16:n_elements(flag)-1] = 1
  return, 1-flag
end

Pro continuum, filename, outfile, bkspace=bkspace
  if n_elements(bkspace) eq 0 then bkspace = 150
  junk = fsc_color(/all, color=ctable)
  h=headfits(filename)
  starname = sxpar(h, 'TARGNAME')
  if size(starname, /tname) eq 'LONG' then starname = sxpar(h, 'OBJECT')
  print, starname
  a=mrdfits(filename, 5, /silent)
  !p.multi=[0,1,2]
  plot, a.wave_opt, a.flux_opt
  mask = maskstar(a.wave_opt, starname) and masktell(a.wave_opt)
  fmask = float(mask)
  w=where(fmask eq 0.0)
  if w[0] ne -1 then fmask[w] = !values.f_nan
  oplot, a.wave_opt, a.flux_opt*fmask*0.9, color=ctable.green
  bset1 = bspline_iterfit(a.wave_opt, a.flux_opt, $
                          invvar=a.ivar_opt*mask, bkspace=bkspace, $
                          yfit=yfit)
  oplot, a.wave_opt, yfit, color=ctable.red
  nflux = a.flux_opt/yfit
  w=where(masktell(a.wave_opt))
  nflux[w] = 1
  plot, a.wave_opt, nflux < 1, yrange=[0, 1.1]
  if n_elements(outfile) ne 0 then begin
     mkhdr, hdr, nflux < 1
     sxaddpar, hdr, 'AIRMASS', sxpar(h, 'AIRMASS')
     sxaddpar, hdr, 'TELESCOP', sxpar(h, 'TELESCOP')
     sxaddpar, hdr, 'INSTRUME', sxpar(h, 'INSTRUME')
     sxaddpar, hdr, 'TELID', sxpar(h, 'TELID')
     sxaddpar, hdr, 'MJD-OBS', sxpar(h, 'MJD-OBS')
     sxaddpar, hdr, 'EQUINOX', sxpar(h, 'EQUINOX')
     sxaddpar, hdr, 'RA', sxpar(h, 'RA')
     sxaddpar, hdr, 'DEC', sxpar(h, 'DEC')
     mwrfits, nflux < 1, outfile, hdr, /create
     mwrfits, a.wave_opt, outfile
  endif
  !p.multi=0
end
