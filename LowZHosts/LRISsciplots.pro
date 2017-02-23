Pro LRISsciplots
  junk = fsc_color(/all, color=ctable, /check_connection)
  bprefixes = ['blue600_4000_d560','blue/blue600_4000_d560','blue600_4000_d460','blue/blue600_4000_d460']
  rprefixes = ['red600_7500','red/red600_7500','red600_5000','red/red600_5000']
  fpostfixes = ['.d560','.d560','.d460','.d460']

  for iprefix = 0, 3 do begin
     bfiles=file_search(bprefixes[iprefix]+'/FinalCR/*_F.fits', count=nfiles)
     for i=0, nfiles-1 do begin
        bflux = x_readspec(bfiles[i], wav=bwave, inflg=2)
        fi=file_info(repstr(bfiles[i], bprefixes[iprefix], rprefixes[iprefix]))
        if ~fi.exists then continue
        rflux = x_readspec(repstr(bfiles[i], bprefixes[iprefix], rprefixes[iprefix]), wav=rwave, inflg=2)
        xrange = minmax([bwave, rwave])
        yrange = [0.0, percentile([bflux, rflux], 0.99)]
        yrange[1] *= 1.1
        name = (reverse(strsplit(bfiles[i], '/', /extract)))[0]
        ps_start, filename=repstr(name, '_F.fits', fpostfixes[iprefix]+'.eps'), /nomatch, /encapsulated
        plot, [0], /nodata, xrange=xrange, yrange=yrange, /xstyle, /ystyle, $
              title=repstr(name, '_F.fits','')
        oplot, bwave, bflux, color=ctable.blue
        oplot, rwave, rflux, color=ctable.red
        ps_end
     endfor
  endfor
end
