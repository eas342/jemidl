Pro plotPSFstars, objs
  resolve_obj, objs[0]
  for iobj=0, n_elements(objs)-1 do begin

     thisDevice=!d.name
     set_plot,'ps'
     device, color=1, encapsulated=0, bits_per_pixel=8
     device, xoffset=0, yoffset=0
     device, xsize=7., ysize=7., /inches
     device, filename='PSFstars/'+objs[iobj]->extract('id')+'PSFstars.ps'

     istars = objs[iobj]->getstars(band='i')
     zstars = objs[iobj]->getstars(band='z')
     ipsf = objs[iobj]->getPSF( band='i', /binfrac, /interp, /modelcen, censize=3 )
     zpsf = objs[iobj]->getPSF( band='z', /binfrac, /interp, /modelcen, censize=3 )
     objs[iobj]->freeimages
     !p.multi=[0,10,10]
     for istar=0, n_elements(istars)-1 do begin
        tvimage, bytscl((bytscl(istars[istar].image/max(istars[istar].image), min=0.00001, max=1))^(1./5)), /nointerp
     endfor
     tvimage, bytscl((bytscl(ipsf, min=0.00001, max=1))^(1./5)), /nointerp
     erase

     !p.multi=[0,10,10]
     for istar=0, n_elements(zstars)-1 do begin
        tvimage, bytscl((bytscl(zstars[istar].image/max(zstars[istar].image), min=0.00001, max=1))^(1./5)), /nointerp
     endfor
     tvimage, bytscl((bytscl(zpsf, min=0.00001, max=1))^(1./5)), /nointerp
     
     device, /close_file
     !p.multi=0
     set_plot, thisDevice
  endfor
end
