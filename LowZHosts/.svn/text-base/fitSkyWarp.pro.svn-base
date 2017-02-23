Function fitSkyWarpFitFunc, p, obj=obj, waveRange=waveRange
  fluxCoeffs = p[0:7] ;; hard-code these for now...
  waveCoeffs = p[8:9]
  param = { waveRange:waveRange, $
            waveCoeffs:waveCoeffs, $
            fluxCoeffs:fluxCoeffs }
  return, obj->modelDeviates( param )
end

Pro fitSkyWarp, fileName, shift
  if n_elements(shift) eq 0 then shift = 0.0d
  junk = fsc_color(/all, color=ctable)

  hdr0=headfits(fileName, /silent)
  targName = strtrim(sxpar(hdr0, 'TARGNAME'), 2)
  expName = strmid(fileName, 7, 9)
  obj=mrdfits(fileName, 5, hdr5, /silent)
  nObj=n_elements(obj)

  if strmid(fileName, 4,1) eq 'r' then waveLimits = [5100.0, 7500.0] $
  else waveLimits = [3150.0, 5600.0]

  ps_start, filename='wavecal/'+expName+'.'+targName+'.skyWarp.ps', $
            /nomatch, xsize=8.0, ysize=10.0, /inches, xoffset=0.25, yoffset=0.5
  openw, lun, 'wavecal/'+expName+'.'+targName+'.skyWarp.dat', /get_lun
  for iObj=0, nObj-1 do begin
     ;; Boxcar Extraction
     sig = (obj[iObj].sky_box-median(obj[iObj].sky_box,3))*sqrt(obj[iobj].sivar_box)
     sigma = biweight_scale(sig,/zero)
     mask = abs(sig) lt 10.0*sigma

     wWave = where(obj[iObj].wave_box gt waveLimits[0] and obj[iObj].wave_box lt waveLimits[1])
     waveRange = minmax(obj[iObj].wave_box[wWave])
     warpObj = obj_new('warpSky', $
                       vacwave=obj[iObj].wave_box[wWave], $
                       skyflux=obj[iObj].sky_box[wWave], $
                       skyMask=mask)
     parinfo = replicate({value:0.d},10)
     parinfo[0].value=1.0
     parinfo[8].value=shift
     boxParms= mpfit('fitSkyWarpFitFunc', $
                     functargs={obj:warpObj, $
                                waveRange:waveRange}, $
                     parinfo=parinfo)
     boxWaveCoeffs = boxParms[8:9]
     boxFluxCoeffs = boxParms[0:7]
     printf, lun, '# BOX '+string(iObj, f='(I02)')
     printf, lun, n_elements(boxWaveCoeffs), boxWaveCoeffs, waveRange, $
             f='(I02,'+string(n_elements(boxWaveCoeffs),f='(I)')+'G20.12,2F12.5)'
     printf, lun, n_elements(boxFluxCoeffs), boxFluxCoeffs, $
             f='(I02,'+string(n_elements(boxFluxCoeffs),f='(I)')+'G20.12)'

     ;; Plot Boxcar
     param = { waveRange:waveRange, $
               waveCoeffs:boxparms[8:9], $
               fluxCoeffs:boxparms[0:7] }
     s = warpObj->warpedSky( param )
     m = warpObj->model()

     xrange = minmax(s.vacwave)
     xrange1 = [xrange[0], xrange[0]+0.1*(xrange[1]-xrange[0])]
     yrange = minmax(m.flux[where(m.vacwave gt xrange1[0] and m.vacwave lt xrange1[1])])
     yrange += 0.3*abs(yrange[1]-yrange[0])*[-1,1]
     datayrange = percentile(s.flux[where(s.vacwave gt xrange1[0] $
                                          and s.vacwave lt xrange1[1])], [0.1,0.99])
     datayrange += abs(datayrange[1]-datayrange[0])*0.3*[-1,1]
     yrange[0] = min([yrange[0], datayrange[0]])
     yrange[1] = max([yrange[1], datayrange[1]])
     !p.multi=[0,1,10]
     plot, [0], /nodata, xrange=xrange1, yrange=yrange, xmargin=[0,0], ymargin=[0,0], $
           title=expName+'.'+targName+' BOX '+string(iObj,f='(I02)'), $
           ytickname=replicate(' ',30), xstyle=1
     oplot, m.vacwave, m.flux, color=ctable.blue
     oplot, s.vacwave, s.flux, ps=10
     for i=1, 9 do begin
        xrange1 = [xrange[0]+0.1*i*(xrange[1]-xrange[0]), $
                   xrange[0]+0.1*(i+1)*(xrange[1]-xrange[0])]
        yrange = minmax(m.flux[where(m.vacwave gt xrange1[0] and m.vacwave lt xrange1[1])])
        yrange += 0.3*abs(yrange[1]-yrange[0])*[-1,1]
        datayrange = percentile(s.flux[where(s.vacwave gt xrange1[0] $
                                             and s.vacwave lt xrange1[1])], [0.1,0.99])
        datayrange += abs(datayrange[1]-datayrange[0])*0.3*[-1,1]
        yrange[0] = min([yrange[0], datayrange[0]])
        yrange[1] = max([yrange[1], datayrange[1]])
        plot, [0], /nodata, xrange=xrange1, yrange=yrange, $
              ymargin=[0,0], xmargin=[0,0], ytickname=replicate(' ',30), xstyle=1
        oplot, m.vacwave, m.flux, color=ctable.blue
        oplot, s.vacwave, s.flux, ps=10
     endfor
     !p.multi=0

     ;; Optimal Extraction
     sig = (obj[iObj].sky_opt-median(obj[iObj].sky_opt,3))*sqrt(obj[iobj].sivar_opt)
     sigma = biweight_scale(sig,/zero)
     mask = abs(sig) lt 10.0*sigma

     wWave = where(obj[iObj].wave_opt gt waveLimits[0] and obj[iObj].wave_opt lt waveLimits[1])
     waveRange = minmax(obj[iObj].wave_opt[wWave])
     warpObj->SetProperty, vacwave=ptr_new(obj[iObj].wave_opt[wWave])
     warpObj->SetProperty, skyflux=ptr_new(obj[iObj].sky_opt[wWave])
     warpObj->SetProperty, skyMask=ptr_new(mask)
     optParms= mpfit('fitSkyWarpFitFunc', $
                     functargs={obj:warpObj, $
                                waveRange:waveRange}, $
                     parinfo=parinfo)
     optWaveCoeffs = optParms[8:9]
     optFluxCoeffs = optParms[0:7]
     printf, lun, '# OPT '+string(iObj, f='(I02)')
     printf, lun, n_elements(optWaveCoeffs), optWaveCoeffs, waveRange, $
             f='(I02,'+string(n_elements(optWaveCoeffs),f='(I)')+'G20.12,2F12.5)'
     printf, lun, n_elements(optFluxCoeffs), optFluxCoeffs, $
             f='(I02,'+string(n_elements(optFluxCoeffs),f='(I)')+'G20.12)'

     ;; Plot Optimal
     param = { waveRange:waveRange, $
               waveCoeffs:optparms[8:9], $
               fluxCoeffs:optparms[0:7] }
     s = warpObj->warpedSky( param )
     m = warpObj->model()

     xrange = minmax(s.vacwave)
     xrange1 = [xrange[0], xrange[0]+0.1*(xrange[1]-xrange[0])]
     yrange = minmax(m.flux[where(m.vacwave gt xrange1[0] and m.vacwave lt xrange1[1])])
     yrange += 0.3*abs(yrange[1]-yrange[0])*[-1,1]
     datayrange = percentile(s.flux[where(s.vacwave gt xrange1[0] $
                                          and s.vacwave lt xrange1[1])], [0.1,0.99])
     datayrange += abs(datayrange[1]-datayrange[0])*0.3*[-1,1]
     yrange[0] = min([yrange[0], datayrange[0]])
     yrange[1] = max([yrange[1], datayrange[1]])
     !p.multi=[0,1,10]
     plot, [0], /nodata, xrange=xrange1, yrange=yrange, xmargin=[0,0], ymargin=[0,0], $
           title=expName+'.'+targName+' OPT '+string(iObj,f='(I02)'), $
           ytickname=replicate(' ',30), xstyle=1
     oplot, m.vacwave, m.flux, color=ctable.blue
     oplot, s.vacwave, s.flux, ps=10
     for i=1, 9 do begin
        xrange1 = [xrange[0]+0.1*i*(xrange[1]-xrange[0]), $
                   xrange[0]+0.1*(i+1)*(xrange[1]-xrange[0])]
        yrange = minmax(m.flux[where(m.vacwave gt xrange1[0] and m.vacwave lt xrange1[1])])
        yrange += 0.3*abs(yrange[1]-yrange[0])*[-1,1]
        datayrange = percentile(s.flux[where(s.vacwave gt xrange1[0] $
                                             and s.vacwave lt xrange1[1])], [0.1,0.99])
        datayrange += abs(datayrange[1]-datayrange[0])*0.3*[-1,1]
        yrange[0] = min([yrange[0], datayrange[0]])
        yrange[1] = max([yrange[1], datayrange[1]])
        plot, [0], /nodata, xrange=xrange1, yrange=yrange, $
              ymargin=[0,0], xmargin=[0,0], ytickname=replicate(' ',30), xstyle=1
        oplot, m.vacwave, m.flux, color=ctable.blue
        oplot, s.vacwave, s.flux, ps=10
     endfor
     !p.multi=0

     ;; Clean up
     obj_destroy, warpObj
  endfor
  ps_end
  close, lun
  free_lun, lun
end
