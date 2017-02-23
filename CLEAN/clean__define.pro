Function CLEAN::Init, $
   image, $
   psf, $
   fwhm, $
   _extra=extra
  
  self.image = ptr_new(image)
  self.PSF = ptr_new(PSF/total(PSF))
  self.workimage = ptr_new(image)
  self.PSFimage = ptr_new(image*0.d)
  self.gaussimage = ptr_new(image*0.d)
  self.resid = ptr_new(image*0.d)
  self.cleanedimage = ptr_new(image*0.d)
  
  nx = (size(psf,/dim))[0]
  ny = (size(psf,/dim))[1]
  self.PSFloc = [(nx-1)/2,(ny-1)/2]
  xs = indgen(nx)
  ys = indgen(ny)
  make_2d, xs, ys
  sig = fwhm/(2.d*sqrt(2.d*alog(2.d)))
  u = ((xs-((nx-1)/2))/sig)^2+((ys-((ny-1)/2))/sig)^2
  
  gauss = exp(-0.5*u)
  gauss /= total(gauss)
  self.gauss = ptr_new(gauss)

  self->SetProperty, _extra=extra
  return, 1
end

Function CLEAN::Clean, maxiter, absfrac=absfrac, absval=absval, verbose=verbose, _extra=extra
  if n_elements(absfrac) eq 0 then absfrac=-1.d
  if n_elements(absval) eq 0 then absval = -10.d
  for iter=0L, maxiter-1 do begin
     if iter mod 100 eq 0 then begin
        if max(*self.workimage)/mean(abs(*self.workimage)) lt absfrac $
        or max(*self.workimage) lt absval then break
        if keyword_set(verbose) then self->print, _extra=extra
     endif
     self->iter1
  endfor
  *self.resid = *self.image - *self.PSFimage
  *self.cleanedimage = *self.gaussimage + *self.resid
  return, (*self.cleanedimage)
end

Function CLEAN::stats
  avgabs = mean(abs(*self.workimage))
  var = (moment(*self.workimage))[1]
  max = max(*self.workimage)
  return, {max:max, avgabs:avgabs, var:var}
end

Pro CLEAN::iter1
  mx=max(*self.workimage,m)  ;; find max pixel
  w=array_indices(*self.workimage,m)  ;; w is location(x/y) of pixel
  ;; add PSF to PSFimage
  *self.PSFimage = AddImage(*self.PSFimage,(*self.PSF)*mx*self.loopgain,w,self.PSFloc)
  ;; add gaussian to gaussimage
  *self.gaussimage = AddImage(*self.gaussimage,(*self.gauss)*mx*self.loopgain,w,self.PSFloc)
  ;; subtract PSF from workimage
  *self.workimage = AddImage(*self.workimage,-1.*(*self.PSF)*mx*self.loopgain,w,self.PSFloc)
end

Pro CLEAN::print, phot=phot
  atv, [[*self.workimage, *self.gaussimage, *self.cleanedimage], $
        [*self.image,*self.PSFimage,*self.resid]]
;  s=size(*self.gaussimage, /dim)
;  atv, rebin([[*self.gaussimage, *self.cleanedimage, *self.cleanedimage*0], $
;              [*self.image,*self.PSFimage,*self.resid]], 2*3*s[0], 2*2*s[1], /sample)
end

Function CLEAN::Phot, center, radius
  return, apphot(*self.gaussimage+*self.workimage, center[0], center[1], radius)
end

Pro CLEAN::CleanUp
  ptr_free, self.image, self.workimage, self.PSFimage, self.gaussimage, self.resid, self.cleanedimage, self.PSF, self.gauss
end

Pro CLEAN__Define, struct
  struct = { CLEAN, $
             Inherits JEMobject, $
             image:        ptr_new(), $
             workimage:    ptr_new(), $
             PSFimage:     ptr_new(), $
             gaussimage:   ptr_new(), $
             resid:        ptr_new(), $
             cleanedimage: ptr_new(), $
             PSF:          ptr_new(), $
             PSFloc:       dblarr(2), $
             gauss:        ptr_new(), $
             fwhm:         0.d, $
             loopgain:     0.d }
end
