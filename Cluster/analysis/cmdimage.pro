Pro cmdimage, obj, filehead
  resolve_obj, obj
  obj->LoadImage
  s=obj->summary()
  z=obj->extract('zcluster')
  colors = obj->color(/xpsf)
  mags=s.zmag_auto

  ebv = obj->extract('ebv')
  A_i = 1.973*ebv
  A_z = 1.472*ebv

  mags -= A_z
  colors -= (A_i - A_z)

  mags += mags*(-0.0422781)+0.635273

  rimage = fltarr(2500,2500)
  gimage = fltarr(2500,2500)
  bimage = fltarr(2500,2500)

  starcheck = s.zflux_radius lt 2.2
  w=where(~starcheck and $
          mags gt 19 and mags lt 25 and $
          colors gt -0.3 and colors lt 1.3 and $
          s.a_image*s.kron_radius lt 200 and $
          s.chi2 ne -1. and s.chi2perdof le 2.)

  for igal=0, n_elements(w)-1 do begin
     counter, igal+1, n_elements(w)
     s1 = s[w[igal]]
     xcenter = fix(2500*(mags[w[igal]]-19)/(25-19))
     ycenter = fix(2500*(colors[w[igal]]-(-0.3))/(1.3-(-0.3)))
     iimage = obj->galimage( w[igal], band='i' )
     zimage = obj->galimage( w[igal], band='z' )
;     mask = zimage*0
     nx = (size(zimage, /dim))[0]
     ny = (size(zimage, /dim))[1]
     dist_ellipse, ellipse, [nx, ny], $
                   s1.xwin_image-s1.xmin-1, s1.ywin_image-s1.ymin-1, $
                   1./s1.boa, s1.pa
;     wmask = where( ellipse lt s1.Re*1.5+10 )
     mask = ellipse lt ((s1.qprad+5) > 5)
;     wmask = where( ellipse lt (s1.qprad+5 > 5) )
;     mask[wmask] = 1
     zimage *= mask
     iimage *= mask
     r = zimage
     g = 0.5*(zimage+iimage)
     b = iimage
     rimage = addimage( rimage, r>0.000, [xcenter,ycenter], [nx/2,ny/2] )
     gimage = addimage( gimage, g>0.000, [xcenter,ycenter], [nx/2,ny/2] )
     bimage = addimage( bimage, b>0.000, [xcenter,ycenter], [nx/2,ny/2] )
  endfor
  obj->FreeImages
  aimage = uindgen(2500,2500)*0L+65535L
  image = uindgen(4,2500,2500)
  image[0,*,*] = 256L*bytscl(rimage, min=0., max=0.09)
  image[1,*,*] = 256L*bytscl(gimage, min=0., max=0.1)
  image[2,*,*] = 256L*bytscl(bimage, min=0., max=0.1)
  image[3,*,*] = aimage
  write_png, filehead+'.png', image
;  mwrfits, rimage, filehead+'R.fits', /create
;  mwrfits, gimage, filehead+'G.fits', /create
;  mwrfits, bimage, filehead+'B.fits', /create
end
