Pro measurePeakCenters, filename, shift
  common measurePeakCentersPersist, skyModel, skyLines

  if n_elements(shift) eq 0 then shift=0.0
  hdr0=headfits(filename, /silent)
  targName = strtrim(sxpar(hdr0, 'TARGNAME'), 2)
  exptime = sxpar(hdr0, 'ELAPTIME')
  if strtrim(sxpar(hdr0, 'INSTRUME'),2) eq 'OSIRIS' then begin
     exptime = sxpar(hdr0, 'EXPTIME')
     targName = strtrim(sxpar(hdr0, 'OBJECT'), 2)
  endif
  targName = lowzhost_nametranslate(targName)
  if strmid(targName,0,2) eq 'HD' then return
  if exptime lt 10 then return
;; from right after 'std-' or 'sci-' to first '.'
  expName = strmid((strsplit(filename, '.', /extract))[0],4)
  obj=mrdfits(filename, 5, hdr5, /silent)
  nObj=n_elements(obj)
  if n_elements(skyModel) eq 0 then begin
     dir=getenv('JEM_REF')+'SPECTRA/SKY/'
     if strtrim(sxpar(hdr0, 'INSTRUME'),2) eq 'OSIRIS' then begin
        skyModel=mrdfits(dir+'sky.fwhm.6.0.fits', 1, /silent)
     endif else begin
        skyModel=mrdfits(dir+'sky.fwhm.3.0.fits', 1, /silent)
     endelse
     readcol, dir+'skylines.dat', skyLines, format='F', /silent
  endif
  lineWindow=10.0 ; Angstroms

  junk=fsc_color(/all, color=ctable, /check_connection)
  if ~(file_info('wavecal')).directory then $
     file_mkdir, 'wavecal'
  ps_start, filename='wavecal/'+targName+'.'+expName+'.peakCenters.ps', /nomatch, /landscape
  openw, lun, 'wavecal/'+targName+'.'+expName+'.peakCenters.dat', /get_lun
  thick=0.5
  for id=0, nObj-1 do begin
;; Boxcar Extraction
     printf, lun, '# BOX ', string(id,f='(I02)')
     waveBox = obj[id].wave_box + shift
     waveBoxRange = minmax(waveBox)

     ;; first refine shift with the 5578 line.
     w=where(waveBox gt 5578.9-15 and waveBox lt 5578.9+15)
     xpeak = find_nminima(-obj[id].sky_box[w], waveBox[w], nfind=1, width=10)
     shift += (5578.9-xpeak)
     waveBox = obj[id].wave_box + shift
     waveBoxRange = minmax(waveBox)

     skywaves = fltarr(n_elements(skyLines))*!values.f_nan
     resids = skywaves
     !p.multi=[0,6,8]
     for iSkyLine=0, n_elements(skyLines)-1 do begin
        ;; determine if data wavelengths cover this skyline
        if (waveBoxRange[0] gt skyLines[iSkyLine]-lineWindow or $
            waveBoxRange[1] lt skyLines[iSkyLine]+lineWindow) then continue
        ;; fit model skyline peak wavelength
        wModel = where(skyModel.vacWave gt skyLines[iSkyLine]-lineWindow/2.0 and $
                       skyModel.vacWave lt skyLines[iSkyLine]+lineWindow/2.0)
        modelFit = mpfitpeak(skyModel.vacWave[wModel], skyModel.flux[wModel], modelFitParam)
        ;; fit data skyline peak wavelength
        wData = where(waveBox gt skyLines[iSkyLine]-lineWindow/2.0 and $
                      waveBox lt skyLines[iSkyLine]+lineWindow/2.0)
        dataFit = mpfitpeak(waveBox[wdata], obj[id].sky_box[wdata], datafitparam)
        skywaves[iSkyLine] = modelfitparam[1]
        resids[iSkyLine] = (datafitparam[1]-shift) - modelfitparam[1]
        if dataFitParam[0] le 0.0 or abs(alog(dataFitParam[2]/modelFitParam[2])) gt alog(2.0) then $
           resids[iSkyLine] = !values.f_nan
        ;; plot results
        xrange = skyLines[iSkyLine]+[-2.5,2.5]*lineWindow
        yrange = minmax(obj[id].sky_box[wData]-dataFitParam[3])
        yrange += abs(yrange[1]-yrange[0])*0.2*[-1,1]
        plot, waveBox, obj[id].sky_box-dataFitParam[3], ps=10, $
              xrange=xrange, yrange=yrange, $
              xstyle=1, ystyle=1, $
              xtickname=replicate(' ',30), $
              ytickname=replicate(' ',30), $
              thick=thick, $
              xmargin=[0,0], ymargin=[0,0]
        x=findgen(200)/199.*(xrange[1]-xrange[0])+xrange[0]
        u=(x-dataFitParam[1])/dataFitParam[2]
        if ~finite(resids[iSkyLine]) then begin
           oplot, x, dataFitParam[0]*exp(-0.5*u^2), thick=thick, color=ctable.red
           oplot, dataFitParam[1]*[1,1], [-1e20, 1e20], thick=thick, color=ctable.red
        endif else begin
           oplot, x, dataFitParam[0]*exp(-0.5*u^2), thick=thick
           oplot, dataFitParam[1]*[1,1], [-1e20, 1e20], thick=thick
        endelse
;        u=(x-modelFitParam[1])/modelFitParam[2]
;        oplot, x, abs(dataFitParam[0])*exp(-0.5*u^2), color=ctable.blue
        oplot, modelFitParam[1]*[1,1], [-1e20, 1e20], color=ctable.blue, thick=thick
        xyouts, xrange[0]+(xrange[1]-xrange[0])*0.1, $
                yrange[0]+(yrange[1]-yrange[0])*0.85, $
                string(modelFitParam[1], f='(F-9.3)'), charsize=0.5
        xyouts, xrange[0]+(xrange[1]-xrange[0])*0.1, $
                yrange[0]+(yrange[1]-yrange[0])*0.75, $
                string(dataFitParam[1]-modelFitParam[1], f='(F-9.3)'), charsize=0.5
;        oplot, skyModel.vacWave, skyModel.flux*dataFitParam[0]/modelFitParam[0], color=ctable.blue
        printf, lun, format='(F-9.1, 3F12.4)', $
                skyLines[iSkyLine], modelFitParam[1], $
                resids[iSkyLine]+modelFitParam[1], resids[iSkyLine]
     endfor
;; Optimal Extraction
     printf, lun, '# OPT ', string(id,f='(I02)')
     waveOpt = obj[id].wave_opt + shift
     waveOptRange = minmax(waveOpt)

     ;; first refine shift with the 5578 line.
     w=where(waveOpt gt 5578.9-15 and waveOpt lt 5578.9+15)
     xpeak = find_nminima(-obj[id].sky_opt[w], waveOpt[w], nfind=1, width=10)
     shift += (5578.9-xpeak)
     waveOpt = obj[id].wave_opt + shift
     waveOptRange = minmax(waveOpt)

     skywaves = fltarr(n_elements(skyLines))*!values.f_nan
     resids = skywaves
     !p.multi=[0,6,8]
     for iSkyLine=0, n_elements(skyLines)-1 do begin
        ;; determine if data wavelengths cover this skyline
        if (waveOptRange[0] gt skyLines[iSkyLine]-lineWindow or $
            waveOptRange[1] lt skyLines[iSkyLine]+lineWindow) then continue
        ;; fit model skyline peak wavelength
        wModel = where(skyModel.vacWave gt skyLines[iSkyLine]-lineWindow/2.0 and $
                       skyModel.vacWave lt skyLines[iSkyLine]+lineWindow/2.0)
        modelFit = mpfitpeak(skyModel.vacWave[wModel], skyModel.flux[wModel], modelFitParam)
        ;; fit data skyline peak wavelength
        wData = where(waveOpt gt skyLines[iSkyLine]-lineWindow/2.0 and $
                      waveOpt lt skyLines[iSkyLine]+lineWindow/2.0)
        dataFit = mpfitpeak(waveOpt[wdata], obj[id].sky_opt[wdata], datafitparam)
        skywaves[iSkyLine] = modelfitparam[1]
        resids[iSkyLine] = (datafitparam[1]-shift) - modelfitparam[1]
        if dataFitParam[0] le 0.0 or abs(alog(dataFitParam[2]/modelFitParam[2])) gt alog(2.0) then $
           resids[iSkyLine] = !values.f_nan
        ;; plot results
        xrange = skyLines[iSkyLine]+[-2.5,2.5]*lineWindow
        yrange = minmax(obj[id].sky_opt[wData]-dataFitParam[3])
        yrange += abs(yrange[1]-yrange[0])*0.2*[-1,1]
        plot, waveOpt, obj[id].sky_opt-dataFitParam[3], ps=10, $
              xrange=xrange, yrange=yrange, $
              xstyle=1, ystyle=1, $
              xtickname=replicate(' ',30), $
              ytickname=replicate(' ',30), $
              thick=thick, $
              xmargin=[0,0], ymargin=[0,0]
        x=findgen(200)/199.*(xrange[1]-xrange[0])+xrange[0]
        u=(x-dataFitParam[1])/dataFitParam[2]
        if ~finite(resids[iSkyLine]) then begin
           oplot, x, dataFitParam[0]*exp(-0.5*u^2), thick=thick, color=ctable.red
           oplot, dataFitParam[1]*[1,1], [-1e20, 1e20], thick=thick, color=ctable.red
        endif else begin
           oplot, x, dataFitParam[0]*exp(-0.5*u^2), thick=thick
           oplot, dataFitParam[1]*[1,1], [-1e20, 1e20], thick=thick
        endelse
;        u=(x-modelFitParam[1])/modelFitParam[2]
;        oplot, x, abs(dataFitParam[0])*exp(-0.5*u^2), color=ctable.blue
        oplot, modelFitParam[1]*[1,1], [-1e20, 1e20], color=ctable.blue, thick=thick
        xyouts, xrange[0]+(xrange[1]-xrange[0])*0.1, $
                yrange[0]+(yrange[1]-yrange[0])*0.85, $
                string(modelFitParam[1], f='(F-9.3)'), charsize=0.5
        xyouts, xrange[0]+(xrange[1]-xrange[0])*0.1, $
                yrange[0]+(yrange[1]-yrange[0])*0.75, $
                string(dataFitParam[1]-modelFitParam[1], f='(F-9.3)'), charsize=0.5
;        oplot, skyModel.vacWave, skyModel.flux*dataFitParam[0]/modelFitParam[0], color=ctable.blue
        printf, lun, format='(F-9.1, 3F12.4)', $
                skyLines[iSkyLine], modelFitParam[1], $
                resids[iSkyLine]+modelFitParam[1], resids[iSkyLine]
     endfor
     !p.multi=0
  endfor
  close, lun
  free_lun, lun
  ps_end
end
