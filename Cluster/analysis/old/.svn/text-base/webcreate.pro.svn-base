Pro webcreate, objs, outdir, _extra=extra
  resolve_obj, objs[0]
  resolve_routine, 'spectrum__define', /com, /no
  resolve_routine, 'dopplershift__define', /com, /no

  readcol, '/home/jmeyers314/supernova.txt', $
           SNname, SNnickname, SNra, SNdec, SNz, $
           format='A,A,X,X,A,A,D,X,X,X,X,X'
  SNgalid = intarr( n_elements(SNname) )
  SNcluster = strarr( n_elements(SNname) )
  specgals = specgals(objs)

  build_starflux, starflux, starwave=starwave, /gal

  ;;;;;;;;;;;;;;;;;;
  ;; load color table

  colors = Obj_New("IDLgrPalette")
  colors->LoadCT, 0
  colors->GetProperty, Red=r, Green=g, Blue=b
  Obj_Destroy, colors
  TVLCT, r, g, b
  red = getcolor('red', !d.table_size-2)
  green = getcolor('green', !d.table_size-3)
  blue = getcolor('blue', !d.table_size-4)

  ;;;;;;;;;;;;;;;;
  ;; process output directory
  
  if strmid(outdir, strlen(outdir)-1,1) ne '/' then outdir += '/'
  print, outdir
  file_mkdir, outdir

  ;;;;;;;;;;;;;;;;;;;
  ;; get cluster names

  clusters=strarr(n_elements(objs))
  for iobj=0, n_elements(objs)-1 do clusters[iobj]=objs[iobj]->Extract('clustername')

  ;;;;;;;;;;;;;;;;;;;
  ;; identify SNe gals
  for iSN=0, n_elements(SNname)-1 do begin
     get_coords, coords, instring=SNra[iSN]+' '+SNdec[iSN]
     alpha=coords[0]*15.d
     delta=coords[1]
     clustername = strmid(SNname[iSN],3,1)
     wcluster = where( clusters eq clustername )
     sexcat = objs[wcluster]->Extract('sexcat')

     dist2 = (sexcat.alphawin_j2000 - alpha)^2*cos(delta*!dpi/180.d)^2 $
             + (sexcat.deltawin_j2000 - delta)^2
     md2 = min( dist2, m)*3600.d^2
     SNgalid[iSN] = m
     SNcluster[iSN] = clustername
  endfor

  ;;;;;;;;;;;;;;;;;;;
  ;; setup for indiv cluster

  currentobj=0
  objs[currentobj]->GetProperty, clustername=currentcluster, zmean=zcluster
  clusterbegin = 0
  clusterend = n_elements( where( specgals.clustername eq currentcluster ) )-1
  file_mkdir, outdir+currentcluster
  objs[currentobj]->LoadiImage
  objs[currentobj]->LoadzImage
  spectra = objs[currentobj]->Extract('spectra')
  window, 10, xsize=800, ysize=350, /pixmap
;  plot, [0], [0], /nodata, xstyle=1, ystyle=1, xrange=[19,25], yrange=[-1,2]
  mags = (objs[currentobj]->Extract('sexcat')).zmag_best
  colors = objs[currentobj]->color( _extra=extra )
;  plotsym, 0, 1.0
;  oplot, mags, colors, ps=8, color=green
;  wspec = where(specgals.clustername eq currentcluster)
;  wwspec = where(is_member((objs[currentobj]->Extract('sexcat')).galid, specgals[wspec].galid))
;  plotsym, 0, 1.0, /fill, color=green
;  oplot, mags[wwspec], colors[wwspec], ps=8
  plot1CMD, objs[currentobj], 1.0, _extra=extra
  

  ;;;;;;;;;;;;;;;;;;;;
  ;; start clusterwide HTML file

  openw, lun, outdir+currentcluster+'/index.html', /get_lun
  printf, lun, '<html><body bgcolor="black"><body text="green">'
  printf, lun, '<title> '+currentcluster+' cluster</title>'

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; loop through spectroscopic gals
  
  for i=0, n_elements(specgals)-1 do begin
     igal=specgals[i]
     print, igal.clustername+string(igal.galid, format='(I04)')

     ;;;;;;;;;;;;;;;;;
     ;; load new cluster if necessary
     if igal.clustername ne currentcluster then begin
        objs[currentobj]->FreeImages
        currentcluster = igal.clustername 
        clusterbegin = i
        clusterend += n_elements( where( specgals.clustername eq currentcluster ) )
        file_mkdir, outdir+currentcluster
        currentobj = (where(clusters eq currentcluster))[0]
        objs[currentobj]->LoadiImage
        objs[currentobj]->LoadzImage
        spectra = objs[currentobj]->Extract('spectra')

        printf, lun, '</html>'
        close, lun
        free_lun, lun
        openw, lun, outdir+currentcluster+'/index.html', /get_lun
        printf, lun, '<html><body bgcolor="black"><body text="green">'
        printf, lun, '<title> '+currentcluster+' cluster</title>'

        wdelete, 10
        window, 10, xsize=800, ysize=350, /pixmap
;        plot, [0], [0], /nodata, xstyle=1, ystyle=1, xrange=[19,25], yrange=[-1,2], $
;              xtitle='z_850', ytitle='i_775 - z_850', title=currentcluster
        mags = (objs[currentobj]->Extract('sexcat')).zmag_best
        colors = objs[currentobj]->color( _extra=extra )
;        plotsym, 0, 1.0
;        oplot, mags, colors, ps=8, color=green
;        wspec = where(specgals.clustername eq currentcluster)
;        wwspec = where(is_member((objs[currentobj]->Extract('sexcat')).galid, specgals[wspec].galid))
;        plotsym, 0, 1.0, /fill, color=green
;        oplot, mags[wwspec], colors[wwspec], ps=8
        plot1CMD, objs[currentobj], 1.0, _extra=extra
     endif

     w = (where( spectra.galid eq igal.galid ))
     if w[0] ne -1 then specs = spectra[w].spectrum $
     else delvarx, specs

     ;;;;;;;;;;;;;;;;;;;;;;;;;
     ; grab ACS image sections

     sectionsize=2*50+1.d
     section = [igal.xwin_image-(sectionsize-1)/2, $
                igal.ywin_image-(sectionsize-1)/2, $
                igal.xwin_image+(sectionsize-1)/2, $
                igal.ywin_image+(sectionsize-1)/2]

     ;;i band
     window, 0, xsize=sectionsize, ysize=sectionsize, /pixmap
     iimagesection = objs[currentobj]->ImageSection(section, band='i')
     tvimage, bytscl(iimagesection, max=max(iimagesection)/2. < 0.1,top=251 ), /nointerp
     for ispec=0, n_elements(specs)-1 do begin
        spec=specs[ispec]
        if obj_valid(spec) then begin
           dx = (spectra[w].ra - igal.alphawin_j2000)*cos(spectra[w].dec*!dpi/180.d)*3600.*20
           dy = (spectra[w].dec - igal.deltawin_j2000)*3600.*20
           plots, 50+dx-5, 50+dy, color=green, /device
           plots, 50+dx+5, 50+dy, color=green, /device, /continue
           plots, 50+dx, 50+dy-5, color=green, /device
           plots, 50+dx, 50+dy+5, color=green, /device, /continue
        endif
     endfor
     ;check for SNe
     w2 = where( SNcluster eq igal.clustername and SNgalid eq igal.galid )
     if w2[0] ne -1 then begin
        get_coords, coords, instring=SNra[w2[0]]+' '+SNdec[w2[0]]
        alpha=coords[0]*15.d
        delta=coords[1]
        snx = (alpha-igal.alphawin_j2000)*cos(delta*!dpi/180.d)*3600.*20
        sny = (delta-igal.deltawin_j2000)*3600.*20
        plots, 50+snx-5, 50+sny, color=red, /device
        plots, 50+snx+5, 50+sny, color=red, /device, /continue
        plots, 50+snx, 50+sny-5, color=red, /device
        plots, 50+snx, 50+sny+5, color=red, /device, /continue
     endif

     void = tvread(filename=outdir+igal.clustername+'/' $
                   +igal.clustername+string(igal.galid,format='(I04)')+'i',/jpeg, /nodialog)
     
     wdelete, 0
     ;;z band
     window, 1, xsize=sectionsize, ysize=sectionsize, /pixmap
     zimagesection = objs[currentobj]->ImageSection(section, band='z')
     tvimage, bytscl(zimagesection, max=max(zimagesection)/2. < 0.1,top=251 ), /nointerp
     for ispec=0, n_elements(specs)-1 do begin
        spec=specs[ispec]
        if obj_valid(spec) then begin
           dx = (spectra[w].ra - igal.alphawin_j2000)*cos(delta*!dpi/180.d)*3600.*20
           dy = (spectra[w].dec - igal.deltawin_j2000)*3600.*20
           plots, 50+dx-5, 50+dy, color=green, /device
           plots, 50+dx+5, 50+dy, color=green, /device, /continue
           plots, 50+dx, 50+dy-5, color=green, /device
           plots, 50+dx, 50+dy+5, color=green, /device, /continue
        endif
     endfor
     if w2[0] ne -1 then begin
        get_coords, coords, instring=SNra[w2[0]]+' '+SNdec[w2[0]]
        alpha=coords[0]*15.d
        delta=coords[1]
        snx = (alpha-igal.alphawin_j2000)*cos(spectra[w].dec*!dpi/180.d)*3600.*20
        sny = (delta-igal.deltawin_j2000)*3600.*20
        plots, 50+snx-5, 50+sny, color=red, /device
        plots, 50+snx+5, 50+sny, color=red, /device, /continue
        plots, 50+snx, 50+sny-5, color=red, /device
        plots, 50+snx, 50+sny+5, color=red, /device, /continue
     endif
     void = tvread(filename=outdir+igal.clustername+'/' $
                   +igal.clustername+string(igal.galid,format='(I04)')+'z',/jpeg, /nodialog)
     wdelete, 1

     ;;;;;;;;;;;;;;;;;;;;;;;;
     ; grab spectrum
     window, 2, xsize=800, ysize=300, /pixmap
     first=1
     delvarx, binbounds
     for ispec=0, n_elements(specs)-1 do begin
        spec=specs[ispec]
        if obj_valid(spec) then begin
           wave=spec->wavelength(frame='obs')
           flux=spec->flux(frame='obs')
           ivar=spec->ivar(frame='obs')
           z=(spec->new_z())->z()
           if first then begin
              jem_bin_spectrum, wave, flux, ivar, dloglam=0.0006, $
                                binflux=binflux, binivar=binivar, binwave=binwave, binbounds=binbounds
              first=0
           endif else begin
              jem_bin_spectrum, wave, flux, ivar, /append, $
                                binflux=binflux, binivar=binivar, binwave=binwave, binbounds=binbounds
           endelse
        endif
     endfor
     if first eq 1 then begin  ;; no valid specs
        xyouts, 400, 150, 'no spec', /device
     endif else begin
        wfinite=where(finite(binivar))
        synflux=dblarr(n_elements(binflux))+1/0
        void = getzchi2(binwave[wfinite], binflux[wfinite], binivar[wfinite], z=z, $
                        starflux=starflux, starwave=starwave, synflux=synflux1 )
        synflux[wfinite]=synflux1
        fsort = sort(binflux)
        nsort = n_elements(fsort)
        yrange = binflux[[fsort[fix(0.03*nsort)], fsort[fix(0.97*nsort)]]]
        yrange += [-0.2, 0.2]*(yrange[1]-yrange[0])
        yrange[0] = yrange[0] < 0.
        xrange = minmax(binwave)
        plot, binwave, binflux, xstyle=9, ystyle=1, xrange=xrange, yrange=yrange, ps=10
        oplot, binwave, 1./sqrt(binivar), color=red, ps=10
        oplot, binwave, synflux, color=blue
        axis, xaxis=1, xstyle=1, xrange=xrange/(1.+z), yrange=yrange
        xyouts, 0.15, 0.7, 'z='+string(z,format='(F6.3)'), /norm
     endelse 
     void = tvread( filename=outdir+igal.clustername+'/' $
                    +igal.clustername+string(igal.galid,format='(I04)')+'spec', /jpeg, /nodialog )
     wdelete,2

     ;;;;;;;;;;;;;;;;;;;;;;;;;
     ; write highlighted CMD

     window, 3, xsize=800, ysize=350, /pixmap
     plot, [0], [0], /nodata, xrange=[19.,25.], yrange=[-0.5,1.5], xstyle=1, ystyle=1
     device, copy=[0,0,800,350,0,0,10]
     w1 = where( (objs[currentobj]->Extract('sexcat')).galid eq igal.galid )
     plotsym, 0, 2.5, /fill, color=red
     oplot, mags[w1], colors[w1], ps=8
     void = tvread( filename=outdir+igal.clustername+'/' $
                    + igal.clustername+string(igal.galid, format='(I04)')+'cmd', /jpeg, /nodialog )
     wdelete, 3

     ;;;;;;;;;;;;;;;;;;;;;;;;;
     ; write galfit image
     infilename = '/home/jmeyers314/scp2/clustergalaxies/' $
                  +igal.clustername $
                  +'/galfit/' $
                  +igal.clustername $
                  +string(igal.galid, format='(I04)') $
                  +'z.fits'
     if (file_info(infilename)).exists ne 0 then begin
        struct = mrdfits(infilename,1,/silent)
;        xsize = (size(struct.image,/dim))[0] > 50
;        ysize = (size(struct.image,/dim))[1] > 50
;        window, 4, xsize=xsize*3, ysize=ysize*2, /pixmap
        msk = struct.intsb*(1-struct.mskim) ge igal.zsky_sigma*1.5
        window, 4, xsize=150, ysize=100, /pixmap
        im = [[struct.errim, struct.intsb*msk, struct.intsb], $
              [struct.image, struct.model, struct.resid]]
        tvimage, bytscl(im, min=-0.02, max=0.16, top=251), /nointerp
        void = tvread( filename=outdir+igal.clustername+'/' $
                       +igal.clustername+string(igal.galid, format='(I04)')+'galfit', $
                       /jpeg, /nodialog )
        wdelete, 4
     endif

     ;;;;;;;;;;;;;;;;;;;;;;;;;
     ; write global HTML

     printf, lun, '   <br>'
     printf, lun, '   <a href="'+igal.clustername+string(igal.galid,format='(I04)')+'.html">' $
             +igal.clustername+string(igal.galid,format='(I04)')+'</a>'

     ;;;;;;;;;;;;;;;;;;;;;;;;;
     ; write indiv HTML
     if i eq clusterbegin then begin 
        pgal = specgals[clusterend]
        ngal = specgals[i+1]
     endif else if i eq clusterend then begin
        pgal = specgals[i-1]
        ngal = specgals[clusterbegin]
     endif else begin
        pgal = specgals[i-1]
        ngal = specgals[i+1]
     endelse

     openw, lun1, outdir+currentcluster+'/'+igal.clustername+string(igal.galid,format='(I04)')+'.html', /get_lun


     printf, lun1, '<html><body bgcolor="black"><body text="green">'
     printf, lun1, '   <title>'+igal.clustername+string(igal.galid,format='(I04)')+'</title>'

     printf, lun1, '   <a href="'+pgal.clustername+string(pgal.galid,format='(I04)')+'.html">' $
             +pgal.clustername+string(pgal.galid,format='(I04)')+'</a>'
     printf, lun1, '   <a href="'+ngal.clustername+string(ngal.galid,format='(I04)')+'.html">' $
             +ngal.clustername+string(ngal.galid,format='(I04)')+'</a>'
     
     printf, lun1, '   <br>'
     printf, lun1, '   <img src="'+igal.clustername+string(igal.galid,format='(I04)')+'i.jpg">'
     printf, lun1, '   <img src="'+igal.clustername+string(igal.galid,format='(I04)')+'z.jpg">'
     printf, lun1, '   <img src="'+igal.clustername+string(igal.galid,format='(I04)')+'galfit.jpg">'
     ; add table
     printf, lun1, '   <table border="1">'
     printf, lun1, '   <tr>'
     printf, lun1, '   <td>redshift</td>'
     printf, lun1, '   <td>'+string(igal.z,format='(F6.3)')+'</td>'
     printf, lun1, '   </tr>'

     printf, lun1, '   <tr>'
     printf, lun1, '   <td>comments</td>'
     printf, lun1, '   <td>'+igal.comment+'</td>'
     printf, lun1, '   </tr>'

     printf, lun1, '   <tr>'
     printf, lun1, '   <td>Asymmetry</td>'
     printf, lun1, '   <td>'+string(igal.asym,format='(F7.4)')+'</td>'
     printf, lun1, '   </tr>'

     printf, lun1, '   <tr>'
     printf, lun1, '   <td>Concentration</td>'
     printf, lun1, '   <td>'+string(igal.conc,format='(F7.4)')+'</td>'
     printf, lun1, '   </tr>'

     printf, lun1, '   <tr>'
     printf, lun1, '   <td>Gini coefficient</td>'
     printf, lun1, '   <td>'+string(igal.gini,format='(F7.4)')+'</td>'
     printf, lun1, '   </tr>'

     clusterra = objs[currentobj]->Extract('ra')
     clusterdec = objs[currentobj]->Extract('dec')
     clusterz = objs[currentobj]->Extract('z')
     radius2 = (clusterra - igal.alphawin_j2000)^2*cos(clusterdec*!dpi/180.d)^2 $
               + (clusterdec - igal.deltawin_j2000)^2
     radius = sqrt(radius2) ;;in degrees
     radius *= !dpi/180.d  ;; in radians
     radius *= lumdist( clusterz, /silent )/(1.+clusterz)^2 ;;in Mpc

     printf, lun1, '   <tr>'
     printf, lun1, '   <td>Cluster radius</td>'
     printf, lun1, '   <td>'+string(radius,format='(F7.4)')+'</td>'
     printf, lun1, '   </tr>'

     printf, lun1, '   </table>'
     ; end table
     printf, lun1, '   <br>'
     printf, lun1, '   <img src="'+igal.clustername+string(igal.galid,format='(I04)')+'cmd.jpg">'
     printf, lun1, '   <br>'
     printf, lun1, '   <img src="'+igal.clustername+string(igal.galid,format='(I04)')+'spec.jpg">'
     printf, lun1, '   <br>'
     printf, lun1, '   <a href="'+pgal.clustername+string(pgal.galid,format='(I04)')+'.html">' $
             +pgal.clustername+string(pgal.galid,format='(I04)')+'</a>'
     printf, lun1, '   <a href="'+ngal.clustername+string(ngal.galid,format='(I04)')+'.html">' $
             +ngal.clustername+string(ngal.galid,format='(I04)')+'</a>'
     
     printf, lun1, '</html>'
     close, lun1
     free_lun, lun1

  endfor
  wdelete, 10
  printf, lun, '</html>'
  close, lun
  free_lun, lun
end
