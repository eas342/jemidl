Pro extract_check, wav, flt, sci, spec, waverange, outfile
  set_plot, 'ps'
  device, filename=outfile, /encapsulated, /color, bits_per_pixel=8
  device, xsize=8, ysize=10, /inches
  device, xoffset=0, yoffset=0
  device, /times, /isolatin1
  loadct, 0
  red = fsc_color('red', 255)
;  wave=interpolate(wav, spec[0].xpos, spec[0].ypos)
  wave=spec[0].wave_opt
  xpixrange = where(wave ge waverange[0] and wave le waverange[1])
  col = mean((spec[0].xpos)[xpixrange])

  pos1=[0.1, 0.77, 0.95, 0.77+0.19]
  pos2=[0.1, 0.53, 0.95, 0.53+0.19]
  pos3=[0.1, 0.29, 0.95, 0.29+0.19]
  pos4=[0.1, 0.05, 0.95, 0.05+0.19]


  plot, [0], /nodata, $
        xrange=waverange, yrange=[-40,40], /xstyle, /ystyle, $
        position=pos1
  fltthumb = transpose(flt[col-40:col+39,min(xpixrange):max(xpixrange)])
  zrange = [0.9, 1.1]
  tvimage, bytscl(fltthumb, min=zrange[0], max=zrange[1], top=254), $
           position=pos1
  oplot, wave, spec[0].xpos-col, color=red

  plot, [0], /nodata, /noerase, xrange=waverange, yrange=[-40,40], /xstyle, /ystyle, $
        position=pos2
  scithumb = transpose(sci[col-40:col+39,min(xpixrange):max(xpixrange)])
  zrange = percentile(scithumb, [0.01, 0.99])
  tvimage, bytscl(scithumb, min=zrange[0], max=zrange[1], top=254), $
           position=pos2
  oplot, wave, spec[0].xpos-col, color=red

  yrange = percentile((spec[0].flux_opt)[xpixrange], [0.01, 0.99])
  yrange += [-1,1]*0.3*(yrange[1]-yrange[0])
  plot, /noerase, spec[0].wave_opt, spec[0].flux_opt, $
        xrange=waverange, yrange=yrange, $
        /xstyle, /ystyle, $
        position=pos3, ps=10

  yrange = percentile((spec[0].flux_box)[xpixrange], [0.01, 0.99])
  yrange += [-1,1]*0.3*(yrange[1]-yrange[0])
  plot, /noerase, spec[0].wave_box, spec[0].flux_box, $
        xrange=waverange, yrange=yrange, $
        /xstyle, /ystyle, $
        position=pos4, ps=10

  device, /close
end

Pro extract_check1, wavfile, fltfile, scifile, waverange, outfile
  wav=mrdfits(wavfile, /silent)
  flt=mrdfits(fltfile, /silent)
  sci=mrdfits(scifile, /silent)
  sky=mrdfits(scifile, 2, /silent)
  spec=mrdfits(scifile, 5, /silent)
  extract_check, wav, flt, sci-sky, spec, waverange, outfile
end
