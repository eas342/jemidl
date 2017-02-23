Pro starresids, objs
  resolve_obj, objs[0]
  loadct, 0
  for iobj=0, n_elements(objs)-1 do begin
     thisDevice=!d.name
     set_plot,'ps'
     device, color=1, encapsulated=0, bits_per_pixel=8
     device, xoffset=0, yoffset=0
     device, xsize=7., ysize=10., /inches
     device, filename='PSFstars/'+objs[iobj]->extract('id')+'resids.ps'

     stars = objs[iobj]->getStars(band='z')
                                ;void =
                                ;objs[iobj]->getPSF(band='z',/binfrac,/interp,/modelcen,censize=3,
                                ;overPSF=overPSF)
     overpsf = objs[iobj]->iterPSF(band='z')
     !p.multi=[0,3,12]
     totals=[0,0]
     for istar=0, n_elements(stars)-1 do begin
        a=fitstar(stars[istar], overPSF)
        totals=[[totals],[a.totals]]
        tvimage, [bytscl(bytscl([a.image,a.model])^(1./5.)),bytscl(a.resid)]
     endfor
     totals=totals[*,1:*]
     objs[iobj]->freeImages
     device,/close_file
     !p.multi=0
     set_plot, thisDevice

     set_plot,'ps'
     device, color=1, encapsulated=0, bits_per_pixel=8
     device, landscape=0
     device, xoffset=0, yoffset=0
     device, xsize=7., ysize=5., /inches
     device, filename='PSFstars/'+objs[iobj]->extract('id')+'flux.ps'
     diff = totals[0,*]/totals[1,*]
     dmag = -2.5*alog10(diff)
     mag = -2.5*alog10(totals[0,*])+24.867
     plot, mag, dmag, xtitle='flux in inner 5x5 pixel region (AB mags)', ytitle='psf mag - ap mag', ps=1
     device,/close_file
     set_plot, thisDevice

  endfor
end

