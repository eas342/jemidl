Pro plotEllipse, x, y, axisRatio, majAxisRadius, theta, _extra=extra
  angles = findgen(100)/99.*2*!pi
  xs = cos(angles)*majAxisRadius
  ys = sin(angles)*majAxisRadius*axisRatio
  xrots = xs*cos(theta) + ys*sin(theta)
  yrots = ys*cos(theta) - xs*sin(theta)
  xrots += x
  yrots += y
  plots, xrots, yrots, _extra=extra
end

Pro apcomp, imfile, specfile
  ;; read in files
  a=mrdfits(imfile, 0, hdr0, /silent)
  sz=size(a, /dim)
  hdr1=headfits(specfile, /silent)
  stat=mrdfits(repstr(specfile, '_F', '_galstruct'), 1, /silent)
  sdss01ap = 1.5*(lumdist(0.1, /silent)/1.1^2) / $
             (lumdist(stat[0].z, /silent)/((1.+stat[0].z)^2))

  ;; extract galaxy center from spec list
  ra=ten(sxpar(hdr1,'RA'))*15.0d
  dec=ten(sxpar(hdr1,'DEC'))

  ;; convert to pixel coordinates
  adxy, hdr0, ra, dec, x, y

  ;; determine platescale of LRIS
  if sxpar(hdr1, 'INSTRUME') eq 'LRISBLUE' then platescale=0.135 $
  else begin
     if strmid(sxpar(hdr1, 'DATE'), 10) gt '2009-07-01' then begin
        platescale=0.135
     endif else begin
        platescale=0.211
     endelse
  endelse

  ;; determine position angle
  PA = (sxpar(hdr1, 'ROTPOSN')+90)*!dpi/180.d
  cornerx=0.5*[-1, -1, 1, 1, -1]
  cornery=sdss01ap*[-1, 1, 1, -1, -1]
  cornerxp = cornerx*cos(pa) + cornery*sin(pa)
  corneryp = cornery*cos(pa) - cornerx*sin(pa)
  xslit = cornerxp/3600.0/0.0001+x
  yslit = corneryp/3600.0/0.0001+y

  loadct, 0
  min=percentile(a, 0.01)
  max=percentile(a, 0.99)
  red=fsc_color('red', 255, /check_connection)
  green=fsc_color('green', 254, /check_connection)
  cyan=fsc_color('cyan', 253, /check_connection)
  plot, [0], /nodata, position=[0.0, 0.0, 1.0, 1.0], $
        xrange=[0,sz[0]-1], yrange=[0, sz[1]-1], /xstyle, /ystyle
  tvimage, bytscl(a, min=min, max=max, top=252)
  plotEllipse, x, y, 1, 1.5/3600.0/0.0001, 0.0, color=red
  plotEllipse, x, y, 1, sdss01ap/3600.0/0.0001, 0.0, color=green
  plots, xslit, yslit, color=cyan
end
