Function SkyIm, ctx, scihdr, fltdir, looksubdirs=looksubdirs
  fltdir=directoryify(fltdir)
  if size( ctx, /n_dim ) eq 2 then nslice=1 $
  else nslice = (size( ctx, /dim ) )[2]
  nx = (size( ctx, /dim ))[0]
  ny = (size( ctx, /dim ))[1]
  sky = dblarr( nx, ny )
  time = dblarr( nx, ny )
  j=1
  ndrizim=sxpar(scihdr, 'NDRIZIM')
  for islice=0, nslice-1 do begin
     for i=0, 31 do begin
        if j gt ndrizim then break
        fltname=sxpar(scihdr, 'D'+string(j,format='(I03)')+'DATA')
        exp_time=sxpar(scihdr, 'D'+string(j,format='(I03)')+'DEXP')
        fltext=strmid(fltname,strpos(fltname, ','))
        fltname=strmid(fltname,0,strpos(fltname,'['))
        case fltext of 
           ',1]' : ext=1
           ',2]' : ext=4
        endcase
        if n_elements(looksubdirs) ne 0 then begin
           file=(file_search(fltdir+looksubdirs+fltname))[0]
        endif else file = fltdir+fltname
        h=headfits(file,ext=ext)
        sky_mean=sxpar(h,'SKY_MEAN')
        print, fltname, ext, sky_mean, sky_mean/exp_time
        j += 1
        mask=(ctx[*,*,islice] and 2L^i) eq 2L^i
        sky += mask*sky_mean
        time += mask*exp_time
     endfor
     if j gt ndrizim then break
  endfor
  sky /= time
  return, sky
end
