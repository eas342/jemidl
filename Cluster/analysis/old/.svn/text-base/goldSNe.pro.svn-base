Pro goldSNe, objs
  resolve_obj,objs[0]
  resolve_routine,'spectrum__define',/com, /no
  resolve_routine,'dopplershift__define',/com, /no
  readcol, 'supernova.txt', name, nickname, ra, dec, zs, format='A,A,X,X,A,A,D,X,X,X,X,X'

  colors = Obj_New("IDLgrPalette")
  colors->LoadCT, 3
  colors->GetProperty, Red=r, Green=g, Blue=b
  Obj_Destroy, colors
  TVLCT, r, g, b
  red = fsc_color('red', !d.table_size-2)
  blue = fsc_color('blue', !d.table_size-3)
;  grey = fsc_color('gray', !d.table_size-4)
;  darkgrey = fsc_color('dark gray', !d.table_size-5)

  thisDevice=!d.name
  set_plot, 'ps'
  device, /encapsulated, color=1, bits_per_pixel=8, landscape=0
  device, xsize=8., ysize=10., /inches
  device, filename='goldSpectra.eps'

  w = where( member( name, ['CL-H-005','CL-K-000','CL-K-018','CL-D-000','CL-R-012','CL-U-004','CL-V-006'] ) )
  name=name[w]
  nickname=nickname[w]
  ra=ra[w]
  dec=dec[w]
  zs=zs[w]
  
  clusternames=strarr(25)
  for iobj=0, n_elements(objs)-1 do clusternames[iobj]=objs[iobj]->Extract('clusterid')
  
  build_starflux, starflux, starwave=starwave, /gal

  !p.multi = [0,1,n_elements(name)]
  for iname=0, n_elements(name)-1 do begin
     get_coords, coords, instring=ra[iname]+' '+dec[iname]
     alpha=coords[0]*15.d
     delta=coords[1]
     clustername = strmid(name[iname],3,1)
     wcluster = where( clusternames eq clustername )
     sexcat = objs[wcluster]->Extract('sexcat')

     dist2 = (sexcat.alphawin_j2000 - alpha)^2*cos(delta*!dpi/180.d)^2 $
             + (sexcat.deltawin_j2000 - delta)^2
     md2 = min( dist2, m)*3600.d^2


     ;;now make spectrum plot
     spectra = objs[wcluster]->Extract('spectra')
     wspec = where(spectra.galid eq m)
     first = 1
;     delvarx, binbounds
;     for ispec=0, n_elements(wspec)-1 do begin
;        spec = spectra[wspec[ispec]].spectrum
;        if not obj_valid(spec) then continue
;        if strpos(spectra[wspec[ispec]].filename, 'Subaru') ne -1 then continue
;        print, spectra[wspec[ispec]].filename
;        wave=spec->wavelength(frame='obs')
;        flux=spec->flux(frame='obs')
;        ivar=spec->ivar(frame='obs')
;        z=(spec->new_z())->z()
;        if first then begin
;            
;           jem_bin_spectrum, wave, flux, ivar, dloglam=0.0006, $
;                             binflux=binflux, binivar=binivar, binwave=binwave, $
;                             binbounds=binbounds
;           first = 0
;        endif else begin
;           jem_bin_spectrum, wave, flux, ivar, /append, $
;                             binflux=binflux, binivar=binivar, binwave=binwave, $
;                             binbounds=binbounds
;        endelse
;     endfor
     wgood = where(obj_valid(spectra[wspec].spectrum))
     newspec = DrizCoaddSpec(spectra[wspec[wgood]].spectrum, dloglam=0.0006, frame='rest')
;     blah = getzchi2(binwave, binflux, binivar, z=z, starflux=starflux, starwave=starwave, synflux=synflux, npoly=2)
     flux=newspec->flux()
     wave=newspec->wavelength()
     blah = getzchi2(wave, flux, newspec->ivar(), z=z, starflux=starflux, starwave=starwave, synflux=synflux, npoly=2)

     fsort = sort(flux)
     nsort = n_elements(fsort)
     yrange = flux[[fsort[fix(0.03*nsort)], fsort[fix(0.97*nsort)]]]
     yrange += [-0.2, 0.2]*(yrange[1]-yrange[0])
     yrange[0] = yrange[0] < 0.
     xrange = minmax(wave)
     plot, wave, flux, xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, ps=10
     oplot, wave, 1./sqrt(newspec->ivar()), color=red, ps=10
     oplot, wave, synflux, color=blue, ps=10
     xyouts, (xrange[1]-xrange[0])*0.1+xrange[0], (yrange[1]-yrange[0])*0.8+yrange[0], strmid(name[iname],3), /data
     print, name[iname],' = ', clustername, '.', strcompress(m,/r),' at z='+strcompress(z,/r)
  endfor
  !p.multi = 0
  device, /close_file
  set_plot, thisDevice
end
