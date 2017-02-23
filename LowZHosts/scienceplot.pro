Pro scienceplot, dir
  f=file_search(dir+'*.fits.gz', count=nfiles)
  targs = strarr(nfiles)
  for i=0, nfiles-1 do begin
     hdr=headfits(f[i], /silent)
     targs[i]=strtrim(sxpar(hdr,'TARGNAME'),2)
  endfor
  uniqtargs = targs[uniq(targs,sort(targs))]
  for i=0, n_elements(uniqtargs)-1 do begin
     counter, i+1, n_elements(uniqtargs)
     delvarx, spectra
     xrange = [!values.f_infinity, -!values.f_infinity]
     yrange = xrange
     for j=0, n_elements(targs)-1 do begin
        if targs[j] ne uniqtargs[i] then continue
        a=mrdfits(f[j],5,/silent)
        wave = a.wave_box
        flux = a.flux_box
        err = 1./sqrt(a.ivar_box)
        flux = medclip(flux, 41, 5)
        flux = gaussfold(wave, flux, 7.0)
        xrange[0] = min([xrange[0],wave])
        xrange[1] = max([xrange[1],wave])
        yrange1 = percentile(flux,[0.01,0.99])
        yrange[0] = yrange[0] < yrange1[0]
        yrange[1] = yrange[1] > yrange1[1]
        if n_elements(spectra) eq 0 then spectra = {wave:wave, flux:flux} $
        else spectra = [spectra,{wave:wave, flux:flux}]
     endfor
     yrange += 0.1*(yrange[1]-yrange[0])*[-1,1]
     junk = fsc_color(/all, color=ctable, /check_connection)
     colors = [ctable.black, ctable.blue, ctable.red, ctable.magenta, ctable.cyan, ctable.yellow, ctable.green, ctable.orange]
     ps_start, filename=uniqtargs[i]+'.eps', /encapsulated, /nomatch, xsize=18, ysize=12, /inches
     plot, [0], /nodata, xrange=xrange, yrange=yrange, /xstyle, /ystyle, title=uniqtargs[i]
     for j=0, n_elements(spectra)-1 do begin
        oplot, spectra[j].wave, spectra[j].flux, color=colors[j]
     endfor
;     multiplot, [1,n_elements(spectra)]
;     for j=0, n_elements(spectra)-1 do begin
;        plot, spectra[j].wave, spectra[j].flux, xrange=xrange, yrange=yrange, /xstyle, /ystyle
;        multiplot
;     endfor
;     multiplot, /reset
     ps_end
  endfor
end
