Function waveticks, axis, index, value
  common tracecheck_common, spec

  return, string(interpol(spec[0].wave_opt, $
                          findgen(n_elements(spec[0].wave_opt)), value), format='(F7.1)')
end

Pro tracecheck, dir
  common tracecheck_common


  loadct, 0
  red=fsc_color('red', 255, /check_connection)
  dir=directoryify(dir)
  print, dir
  f1=file_search(dir+'Science_std/*.fits.gz', count=nf1)
  f2=file_search(dir+'Science_box/*.fits.gz', count=nf2)
  if nf1 eq 0 then begin
     if nf2 eq 0 then return
     f=f2
  endif else begin
     if nf2 gt 0 then f=[f1,f2] $
     else f=f1
  endelse
  split = strsplit(dir, '/', /extract)
  split = split[n_elements(split)-2:n_elements(split)-1]
  filename = 'Trace_'+strjoin(split)+'.ps'

  ps_start, filename=filename, xsize=8, ysize=10, $
            xoffset=0, yoffset=0, charsize=0.9, /nomatch
  for i=0, n_elements(f)-1 do begin
     counter, i+1, n_elements(f)
     if strpos(f[i],'newwave') ne -1 then continue
     h=headfits(f[i], /silent)
     targName=sxpar(h, 'TARGNAME')
     if strmid(targName,0,2) eq 'HD' then continue
     name=(reverse(strsplit(f[i], '/', /extract)))[0]
     name=repstr(name, '.fits.gz', '')

     sci=mrdfits(f[i], 0, /silent)
     sky=mrdfits(f[i], 2, /silent)
     spec=mrdfits(f[i], 5, /silent)
     sz=size(sci, /dim)
     im = rebin(sci-sky, sz[0], sz[1]/8)
     wx = median(spec[0].xpos)
     wid = 200
     im = im[wx-wid/2:wx+wid/2,*]
     erase
     tvimage, bytscl(im, min=percentile(im, 0.005), max=percentile(im, 0.995), top=254), $
              position=[0.1, 0.1, 0.9, 0.9]
     xyouts, 0.5, 0.95, targName+' '+name, align=0.5, /normal
     plot, /noerase, spec[0].xpos+1, spec[0].ypos, $
           xrange=[wx-wid/2, wx+wid/2+1], yrange=[0, sz[1]], $
           /xstyle, ystyle=9, color=red, position=[0.1, 0.1, 0.9, 0.9]
     axis, yaxis=1, ytickformat='waveticks', color=red, ystyle=1
  endfor
  ps_end
end

Pro dotracecheck
  dir='~/scp4/LowZHosts/JMreduction/'
  sdirs=file_search(dir+'*', count=nsdirs)
  ssdirs=['red600_7500','red600_5000','blue600_4000_d460','blue600_4000_d560']
  nssdirs=n_elements(ssdirs)
  for i=0, nsdirs-1 do begin
     for j=0, nssdirs-1 do begin
        tracecheck, sdirs[i]+'/'+ssdirs[j]
     endfor
  endfor
end
