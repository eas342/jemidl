Pro sky5575plots, dir
  junk = fsc_color(/all, color=ctable, /check_connection)
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
  filename = 'Sky5575_'+strjoin(split)+'.ps'

  ps_start, filename=filename, xsize=8, ysize=10, $
            xoffset=0, yoffset=0, charsize=0.9, /nomatch
  !p.multi=[0,1,8]
  j=0
  for i=0, n_elements(f)-1 do begin
     counter, i+1, n_elements(f)
     if strpos(f[i],'newwave') ne -1 then continue
     name=(reverse(strsplit(f[i], '/', /extract)))[0]
     name=repstr(name, '.fits.gz', '')
     a=mrdfits(f[i], 5, /silent)
     h=headfits(f[i], /silent)
     targName=sxpar(h, 'TARGNAME')
     w=where(a.wave_opt ge 5450 and a.wave_opt le 5700)
;;     yrange = [0.0, percentile(a.sky_opt[w], 0.99)*1.1]
     yrange = minmax(a.sky_opt[w])
     yrange += [-1,1]*0.1*(yrange[1]-yrange[0])
     plot, a.wave_opt, a.sky_opt, xrange=[5450, 5700], yrange=yrange, $
           yticks=1, yminor=1, ytickname=replicate(' ',30), title=name+' '+targName, $
           /xstyle, /ystyle
     newwavefile = repstr(f[i],'.fits.gz','.newwave.fits.gz')
     if (file_info(newwavefile)).exists then begin
        b=mrdfits(newwavefile,5,/silent)
        oplot, b.wave_opt, b.sky_opt, color=ctable.blue
     endif
     if j mod 8 eq 7 then erase
     j += 1
  endfor
  ps_end
  !p.multi=0
end

Pro doSky5575plots
  dir='~/scp4/LowZHosts/JMreduction/'
  sdirs=file_search(dir+'*', count=nsdirs)
  ssdirs=['red600_7500','red600_5000']
  nssdirs=n_elements(ssdirs)
  for i=0, nsdirs-1 do begin
     for j=0, nssdirs-1 do begin
        sky5575plots, sdirs[i]+'/'+ssdirs[j]
     endfor
  endfor
end
