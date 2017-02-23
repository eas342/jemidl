Pro cmdmovie, obj
  resolve_obj, obj
  obj->LoadiImage
  obj->LoadzImage
  s=obj->summary()
  colors = obj->color(/clncat)
  mags=s.zmag_auto

  ebv = obj->extract('ebv')
  A_i = 1.973*ebv
  A_z = 1.472*ebv

  mags -= A_z
  colors -= (A_i - A_z)

  mags += mags*(-0.0422781)+0.635273

  starcheck = s.zflux_radius lt 2.2
  w=where(~starcheck and $
          mags gt 19 and mags lt 25 and $
          colors gt -0.3 and colors lt 1.3 and $
          s.a_image*s.kron_radius lt 200 and $
          s.chi2 ne -1. and s.chi2perdof le 2.)

;  size = [2500,2500]

  size = size(obj->extract('zimage'),/dim)
  xsize = size[0]/5
  ysize = size[1]/5
;  xsize = 1920
;  ysize = 1080
  xrange = minmax(s.xwin_image)
  yrange = minmax(s.ywin_image)

  ;;load thumbnails...
  thumb = { xstart:0L, $
            ystart:0L, $
            xend:0L, $
            yend:0, $
            nx:0L, $
            ny:0L, $
            r:ptr_new(), $
            g:ptr_new(), $
            b:ptr_new() }
  thumb = replicate(thumb, n_elements(w))
  for igal=0, n_elements(w)-1 do begin
     counter, igal+1, n_elements(w), 'preparing thumbnail: '
     s1 = s[w[igal]]
     thumb[igal].xend = fix(xsize*(mags[w[igal]]-19)/(25-19))
     thumb[igal].yend = fix(ysize*(colors[w[igal]]-(-0.5))/(1.5-(-0.5)))
     thumb[igal].xstart = fix(xsize*(s1.xwin_image-xrange[0])/(xrange[1]-xrange[0]))
     thumb[igal].ystart = fix(ysize*(s1.ywin_image-yrange[0])/(yrange[1]-yrange[0]))
     iimage = obj->galimage( w[igal], band='i')
     zimage = obj->galimage( w[igal], band='z')
     mask = zimage*0
     size=size(iimage,/dim)
     if size[0] lt 10 or size[1] lt 10 then continue
     dist_ellipse, ellipse, size, $
                   s1.xwin_image-s1.xmin-1, s1.ywin_image-s1.ymin-1, $
                   1./s1.boa_guess, s1.pa_guess
     wmask = where( ellipse lt s1.r_e*1.5+10 )
     mask[wmask] = 1
     zimage *= mask
     iimage *= mask
     iimage = undersampleimage(iimage, 5)
     zimage = undersampleimage(zimage, 5)
     thumb[igal].nx = (size(iimage,/dim))[0]
     thumb[igal].ny = (size(iimage,/dim))[1]
     thumb[igal].r = ptr_new(zimage)
     thumb[igal].g = ptr_new(0.5*(iimage+zimage))
     thumb[igal].b = ptr_new(iimage)
  endfor

  ;;iterate frames...
  nframes = 20.
  for iframe=0, nframes-1 do begin
     rimage = fltarr(xsize,ysize)
     gimage = fltarr(xsize,ysize)
     bimage = fltarr(xsize,ysize)
     for igal=0, n_elements(w)-1 do begin
        counter, igal+1, n_elements(w), 'Frame '+str(iframe)+' galaxy: '
        i = iframe/(nframes-1)
        x = fix(thumb[igal].xstart*(1-i) + i*thumb[igal].xend)
        y = fix(thumb[igal].ystart*(1-i) + i*thumb[igal].yend)
        nx = thumb[igal].nx
        ny = thumb[igal].ny
        if ~ptr_valid(thumb[igal].r) then continue
        rimage = addimage( rimage, *(thumb[igal].r)>0.001, [x,y], [nx/2,ny/2] )
        gimage = addimage( gimage, *(thumb[igal].g)>0.001, [x,y], [nx/2,ny/2] )
        bimage = addimage( bimage, *(thumb[igal].b)>0.001, [x,y], [nx/2,ny/2] )
     endfor
     print
     image = uindgen(4,xsize,ysize)
     aimage = uindgen(xsize,ysize)*0L+65535L
     image[0,*,*] = 256L*bytscl(rimage, min=0., max=0.09)
     image[1,*,*] = 256L*bytscl(gimage, min=0., max=0.1)
     image[2,*,*] = 256L*bytscl(bimage, min=0., max=0.1)
     image[3,*,*] = aimage
     write_png, obj->extract('clusterid')+string(iframe,format='(I02)')+'.png', image
  endfor
end
