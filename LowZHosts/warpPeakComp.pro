Pro warpPeakComp
  fitPeaksFiles = file_search('*.fitpeaks.dat')
  junk = fsc_color(/all, color=ctable)
  currentTargName = ' '
  ps_start, filename='warpPeakComp.ps', /nomatch
  for i=0, n_elements(fitPeaksFiles)-1 do begin
     targName = strtrim((strsplit(fitPeaksFiles[i],'.',/extract))[1],2)

     ;; if first exposure of given target, then first plot warp data
     if targName ne currentTargName then begin
        fitSkyFiles = file_search('*'+targName+'*fitsky.dat')
        for j=0, n_elements(fitSkyFiles)-1 do begin
           readcol, fitSkyfiles[j], fluxCoeff0, fluxCoeff1, fluxCoeff2, fluxCoeff3, fluxCoeff4, fluxCoeff5, $
                    waveCoeff0, waveCoeff1, waveRange0, waveRange1, f='F,F,F,F,F,F,F,F,F,F', /silent
           for iLine=0, n_elements(waveCoeff0)-1 do begin
              waveCoeffs = [waveCoeff0[iLine], waveCoeff1[iLine]]
              x = [5000.0, 8000.0]
              xnorm = (x-waveRange0[iLine])/(waveRange1[iLine]-waveRange0[iLine])
              waveResid = waveCoeffs##fchebyshev(xnorm, 2)
              if targName ne currentTargName then begin
                 plot, x, waveResid, xrange=x, yrange=[-1.0, 1.0], /xstyle, /ystyle, $
                       title=targName
              endif else begin
                 oplot, x, waveResid
              endelse
              currentTargName = targName
           endfor
        endfor
     endif
     readcol, fitPeaksFiles[i], waveCoeff0, waveCoeff1, waveRange0, waveRange1, format='X,F,F,F,F,X', /silent
     for iLine=0, n_elements(waveCoeff0)-1 do begin
        waveCoeffs = [waveCoeff0[iLine], waveCoeff1[iLine]]
        x = [5000.0, 8000.0]
        xnorm = (x-waveRange0[iLine])/(waveRange1[iLine]-waveRange0[iLine])
        waveResid = waveCoeffs##fchebyshev(xnorm, 2)
        oplot, x, waveResid, color=ctable.red
     endfor
     currentTargName = targName
  endfor
  ps_end
end
