Pro LickCheck, file
  junk = fsc_color(/all, color=ctable, /check_connection)
  Hdelta = [4041.60,  4079.75,  4083.50,  4122.25,  4128.50,  4161.00]
  Hgamma = [4283.50,  4319.75,  4319.75,  4363.50,  4367.25,  4419.75]
  Hbeta  = [4827.875, 4847.875, 4847.875, 4876.625, 4876.625, 4891.625]
  Mgb    = [5142.625, 5161.375, 5160.125, 5192.625, 5191.375, 5206.375]
  Fe5270 = [5233.15,  5248.15,  5245.65,  5285.65,  5285.65,  5318.15]
  Fe5335 = [5304.625, 5315.875, 5312.125, 5352.125, 5353.375, 5363.375]

  LickBins = [[Hdelta],[Hgamma],[Hbeta],[Mgb],[Fe5270],[Fe5335]]
  LickNames = ['H\delta_A', 'H\gamma_A', 'H\beta_A', 'Mg_b', 'Fe5270', 'Fe5335']

  check = mrdfits(file, 1, /silent)
  sz=size(check.indiv,/dim)
  nobj=sz[0]
  npix=sz[1]
  name=repstr(file,'_check.fits','')
  galstruct = mrdfits(name+'_galstruct.fits',1)
  wave = check.wave/(1.0+galstruct[0].z)
  ps_start, filename=repstr(file,'_check.fits','_Lick.ps'), /nomatch, xsize=8, ysize=10, /inches, xoffset=0, yoffset=0
  ;; full overplot
  !p.multi=[0,2,3]
  for iLick=0, 5 do begin
     xrange = [LickBins[0, iLick], LickBins[5, iLick]]
     w=where(wave ge xrange[0] and wave le xrange[1])
     if w[0] eq -1 then $
        yrange = [0,1] $
     else begin
        w1 = where(~check.mask $
                   and check.indiv ne 0 $
                   and transpose(rebin(wave,npix,nobj)) ge xrange[0] $
                   and transpose(rebin(wave,npix,nobj)) le xrange[1])
        yrange = percentile(check.indiv[w1],[0.04,0.96])
     endelse
     xrange += 0.1*[-1,1]*(xrange[1]-xrange[0])
     yrange += 0.3*[-1,1]*(yrange[1]-yrange[0])
     plot, [0], /nodata, xrange=xrange, yrange=yrange, $
           /xstyle, /ystyle
     axis, xaxis=1, xrange=xrange*(1.0+galstruct[0].z)
     polyfill, [LickBins[0,iLick], LickBins[0,iLick], LickBins[1,iLick], LickBins[1,iLick], LickBins[0,iLick]], $
               [yrange[0], yrange[1], yrange[1], yrange[0], yrange[0]], $
               color=ctable.skyblue
     polyfill, [LickBins[4,iLick], LickBins[4,iLick], LickBins[5,iLick], LickBins[5,iLick], LickBins[4,iLick]], $
               [yrange[0], yrange[1], yrange[1], yrange[0], yrange[0]], $
               color=ctable.skyblue
     polyfill, [LickBins[2,iLick], LickBins[2,iLick], LickBins[3,iLick], LickBins[3,iLick], LickBins[2,iLick]], $
               [yrange[0], yrange[1], yrange[1], yrange[0], yrange[0]], $
               color=ctable.papaya
     xyouts, xrange[0]+0.1*(xrange[1]-xrange[0]), yrange[0]+0.9*(yrange[1]-yrange[0]), textoidl(LickNames[iLick])
     if w[0] eq -1 then continue
     oplot, wave, check.coadd, thick=4, ps=10
     w1 = where(check.indiv eq 0.0) ;; don't plot these points...
     if w1[0] ne -1 then check.indiv[w1] = !values.f_nan
     for iobj=0, nobj-1 do begin
        oplot, wave, check.indiv[iobj,*], thick=0.5, ps=10, color=ctable.darkgreen
        wbad = where(check.mask[iobj,*], nbad)
        if nbad gt 0 then $
           oplot, wave[wbad], check.indiv[iobj,wbad], thick=0.5, ps=4, color=ctable.red
     endfor
  endfor
  erase
  !p.multi=0
  ;; indiv overplots
  for iobj=0, nobj-1 do begin
     !p.multi=[0,2,3]
     for iLick=0, 5 do begin
        xrange = [LickBins[0, iLick], LickBins[5, iLick]]
        w=where(wave ge xrange[0] and wave le xrange[1])
        if w[0] eq -1 then $
           yrange = [0,1] $
     else begin
        w1 = where(~check.mask $
                   and check.indiv ne 0 $
                   and transpose(rebin(wave,npix,nobj)) ge xrange[0] $
                   and transpose(rebin(wave,npix,nobj)) le xrange[1])
        yrange = percentile(check.indiv[w1],[0.04,0.96])
     endelse
        xrange += 0.1*[-1,1]*(xrange[1]-xrange[0])
        yrange += 0.3*[-1,1]*(yrange[1]-yrange[0])
        plot, [0], /nodata, xrange=xrange, yrange=yrange, $
              /xstyle, /ystyle
        axis, xaxis=1, xrange=xrange*(1.0+galstruct[0].z)
        polyfill, [LickBins[0,iLick], LickBins[0,iLick], LickBins[1,iLick], LickBins[1,iLick], LickBins[0,iLick]], $
                  [yrange[0], yrange[1], yrange[1], yrange[0], yrange[0]], $
                  color=ctable.skyblue
        polyfill, [LickBins[4,iLick], LickBins[4,iLick], LickBins[5,iLick], LickBins[5,iLick], LickBins[4,iLick]], $
                  [yrange[0], yrange[1], yrange[1], yrange[0], yrange[0]], $
                  color=ctable.skyblue
        polyfill, [LickBins[2,iLick], LickBins[2,iLick], LickBins[3,iLick], LickBins[3,iLick], LickBins[2,iLick]], $
                  [yrange[0], yrange[1], yrange[1], yrange[0], yrange[0]], $
                  color=ctable.papaya
        xyouts, xrange[0]+0.1*(xrange[1]-xrange[0]), yrange[0]+0.9*(yrange[1]-yrange[0]), textoidl(LickNames[iLick])
        if w[0] eq -1 then continue
        oplot, wave, check.coadd, thick=4, ps=10
        oplot, wave, check.indiv[iobj,*], thick=0.5, ps=10, color=ctable.darkgreen
        wbad = where(check.mask[iobj,*], nbad)
        if nbad gt 0 then $
           oplot, wave[wbad], check.indiv[iobj,wbad], thick=0.5, ps=4, color=ctable.red
     endfor
     !p.multi=0
  endfor
  ps_end
end

Pro doLickCheck
  ;; current directory should be scp4/LowZHosts/JMreduction/
  f=file_search('./*/*/Final/*_check.fits')
  for i=0, n_elements(f)-1 do begin
     lickCheck, f[i]
  endfor
end
