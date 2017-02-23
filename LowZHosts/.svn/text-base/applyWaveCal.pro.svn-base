Pro applyWaveCal, sciFileName, calFileName
  cal = readWaveCal(calFileName)
  data = mrdfits(sciFileName, 5, hdr5, /silent)
  hdr0 = headfits(sciFileName, ext=0, /silent)

  newdata = replicate({old_wave_box:dblarr(n_elements(data[0].wave_box)), $
                       old_wave_opt:dblarr(n_elements(data[0].wave_opt))}, $
                      n_elements(data))
  data = struct_addtags(data, newdata)

  ra=ten(sxpar(hdr0,'RA'))*15.0
  dec=ten(sxpar(hdr0,'DEC'))
  jd=double(sxpar(hdr0,'MJD-OBS'))+2400000.5d
  equinox=sxpar(hdr0,'EQUINOX')
  ut=sxpar(hdr0,'UT')
  if ut eq 0 then ut=sxpar(hdr0,'UTC')
  dv=(-1.0)*x_keckhelio(ra, dec, equinox, jd=jd, obs='Keck')
  lambda_corr = sqrt((1.d + dv/299792.458d)/(1.d - dv/299792.458d))
  print, 'UTC: '+string(ut)
  print, 'Heliocentric Correction: '+string(dv)
  print, 'Lambda factor: '+string(lambda_corr)
  print

  for iExtract=0, n_elements(data)-1 do begin
     ;; boxcar
     w=where(cal.isbox and cal.iObj eq iExtract)
     waveCoeffs = *(cal[w].waveCoeffs)
     waveRange = cal[w].waveRange
     waveUnCal = data[iExtract].wave_box
     waveNorm = (waveUnCal - waveRange[0])/(waveRange[1]-waveRange[0])
     waveCal = waveUnCal - waveCoeffs##fchebyshev(waveNorm, n_elements(waveCoeffs))
     waveHelio = waveCal*lambda_corr
     data[iExtract].old_wave_box = data[iExtract].wave_box
     data[iExtract].wave_box = waveHelio
     ;; optimal
     w=where(cal.isopt and cal.iObj eq iExtract)
     waveCoeffs = *(cal[w].waveCoeffs)
     waveRange = cal[w].waveRange
     waveUnCal = data[iExtract].wave_opt
     waveNorm = (waveUnCal - waveRange[0])/(waveRange[1]-waveRange[0])
     waveCal = waveUnCal - waveCoeffs##fchebyshev(waveNorm, n_elements(waveCoeffs))
     waveHelio = waveCal*lambda_corr
     data[iExtract].old_wave_opt = data[iExtract].wave_opt
     data[iExtract].wave_opt = waveHelio
  endfor
  newFileName = repstr(sciFileName, '.fits.gz', '.newwave.fits')
  sxdelpar, hdr0, 'NAXIS2'
  mwrfits, [0], newFileName, hdr0, /create
  mwrfits, [0], newFileName
  mwrfits, [0], newFileName
  mwrfits, [0], newFileName
  mwrfits, [0], newFileName
  mwrfits, data, newFileName
  spawn, 'gzip -f '+newFileName
end
