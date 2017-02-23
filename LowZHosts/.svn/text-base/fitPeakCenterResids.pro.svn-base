Pro fitPeakCenterResids, fileName
  junk = fsc_color(/all, color=ctable, /check_connection)
  fileSplit = strsplit(fileName, '/', /extract)
  fileSplit = (reverse(fileSplit))[0]
  fileSplit = strsplit(fileSplit, '.', /extract)
  expName = fileSplit[0]
  targName = fileSplit[1]
  readcol, fileName, modelWave, dataWave, format='X,F,F', /silent
  wSec = where(modelWave-shift(modelWave,1) lt 0)
  wSec = [wSec, n_elements(modelWave)]

  datFileName = repstr(fileName, 'peakCenters.dat', 'peakCenterResidFit.dat')
  psFileName = repstr(fileName, 'peakCenters.dat', 'peakCenterResidFit.ps')

  angstrom=string(197B)
  openw, lun, datFileName, /get_lun
  ps_start, filename=psFileName, /nomatch, /landscape
  for iSec=0, n_elements(wSec)-2 do begin
     iExtract = iSec/2
     if iSec mod 2 eq 0 then begin
        printf, lun, '# BOX ', string(iExtract,f='(I02)')
        box=1
        opt=0
     endif else begin
        printf, lun, '# OPT ', string(iExtract,f='(I02)')
        box=0
        opt=1
     endelse

     dataWaveSec = dataWave[wSec[iSec]:wSec[iSec+1]-1]
     modelWaveSec = modelWave[wSec[iSec]:wSec[iSec+1]-1]
     w=where(finite(dataWaveSec) and finite(modelWaveSec))
     xwave=dataWaveSec[w]
     x = (xwave-min(xwave))/(max(xwave)-min(xwave))
     y = dataWaveSec[w] - modelWaveSec[w]

     ;; fit with chebyshev polynomial
     order = 1
     yfit = y ;; initialize yfit
     sigma = 1 ;; initialize sigma
     for ifit=0, 2 do begin
        inliers = abs(y-yfit)/sigma lt 3.0 and xwave lt 7500
        wIn = where(inliers, complement=wOut)
        fitparams = svdfit(x[wIn], y[wIn], order+1, function_name='fchebyshev')
        yfit = fitparams##fchebyshev(x, order+1)
        sigma = biweight_scale((y-yfit)[wIn], /zero)
     endfor

     ;; plotable solution
     xplot = findgen(300)/299.0
     yplot = fitparams##fchebyshev(xplot, order+1)
     xplot = xplot*(max(xwave)-min(xwave))+min(xwave)

     ;; plot!
     !p.multi=[0,1,2]
     xrange = minmax(modelWaveSec)
     xrange += abs(xrange[1]-xrange[0])*0.1*[-1,1]
     yrange = minmax(y[wIn])
     yrange += abs(yrange[1]-yrange[0])*0.5*[-1,1]
     plot, xwave[wIn], y[wIn], ps=1, $
           xrange=xrange, yrange=yrange, $
           xstyle=1, ystyle=1, $
           xtitle=textoidl('\lambda_{observed} ('+angstrom+')'), $
           ytitle=textoidl('\lambda_{observed}-\lambda_{model} ('+angstrom+')')
     if wOut[0] ne -1 then $
        oplot, xwave[wOut], y[wOut], color=ctable.red, ps=2
     oplot, xplot, yplot
     yrange = sigma*[-9.0, 9.0]
     plot, xwave[w], y[w]-yfit[w], ps=1, $
           xrange=xrange, yrange=yrange, $
           xstyle=1, ystyle=1, $
           xtitle=textoidl('\lambda_{observed} ('+angstrom+')'), $
           ytitle=textoidl('Residual')
     if wOut[0] ne -1 then $
        oplot, xwave[wOut], y[wOut]-yfit[wOut], color=ctable.red, ps=2
     oplot, xrange, xrange*0
     xyouts, 0.18, 0.91, targName, /normal, charsize=0.9
     xyouts, 0.18, 0.89, expName, /normal, charsize=0.9
     if box then $
        xyouts, 0.18, 0.87, 'Boxcar', /normal, charsize=0.9 $
     else $
        xyouts, 0.18, 0.87, 'Optimal', /normal, charsize=0.9
     xyouts, 0.18, 0.85, 'Extraction '+strtrim(iExtract,2), /normal, charsize=0.9
     xyouts, 0.18, 0.40, 'rms residual: '+string(sigma,f='(F8.3)'), /normal, charsize=0.9
     printf, lun, format='(I02,'+string(order+1,f='(I)')+'G20.12,2F12.5,F10.5)', $
             order+1, fitparams, min(xwave), max(xwave), sigma
     printf, lun, format='(I02)', 0
  endfor
  ps_end
  close, lun
  free_lun, lun
end
