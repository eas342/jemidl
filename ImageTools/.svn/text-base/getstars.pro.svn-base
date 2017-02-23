;+
; NAME:
;	GETSTARS
;
; PURPOSE:
;       Generate an array of star structures for use with PSF construction.
; CALLING SEQUENCE:
;	Stars = GETSTARS(coords, image, errim=errim)
; INPUTS:
;       coords: a [N,2] array of pixel coordinates of stars
;       (SExtractor coordinate convention)
;       image: science image
;      
; KEYWORD PARAMETERS:
;       errim: optional error image
;       starsize: size of thumbnail to be extracted for each star
;       innersize: size of thumbnail used to model star
;
; OUTPUTS:
;       Stars: an array of structures with the following fields:
;              image: a STARSIZExSTARSIZE image of the star
;              errim: a STARSIZExSTARSIZE error image of the star
;              model: the model fit of the star
;              modelparams: fit parameters for a gaussian (from MPFIT2DPEAK)
;                    0 : offset
;                    1 : amplitude
;                    2 : x-sigma
;                    3 : y-sigma
;                    4 : x-center
;                    5 : y-center
; EXAMPLE:
;
; NOTES:
;
; MODIFICATION HISTORY:
;       JEM, Jan, 2009. Written
;       JEM, Jun, 2010. Rewritten to be more independent of data type
;-

Function GetStars, $
   coords, $
   image, $
   errim=errim, $
   background=background, $
   starsize=starsize, $
   innersize=innersize, $
   names=names

  if n_elements(starsize) eq 0 then starsize=41
  if n_elements(innersize) eq 0 then innersize=9
  halfsize=(starsize-1)/2

  onestar = { name:'', $
              coords: fltarr(2), $
              image:dblarr(starsize,starsize), $
              errim:dblarr(starsize,starsize), $
              model:dblarr(innersize,innersize), $
              modelparams: dblarr(7) }
  sz = size(coords, /dim)
  nstars = sz[0]
  if n_elements(names) eq 0 then names = 'star'+str(lindgen(nstars))
  if n_elements(background) eq 0 then background=fltarr(nstars)
  for i=0, nstars-1 do begin
     onecoord = coords[i,*]
     star = onestar
     star.coords=onecoord
     star.name=names[i]
     p = round([onecoord[0]-1, $ ;;need to take sex coords -> IDL coords
                onecoord[1]-1, $
                onecoord[0]-1, $
                onecoord[1]-1])+[-1,-1,1,1]*halfsize
     star.image = imagesection(image, p)
     if finite(background[i]) then star.image -= background[i]
     if n_elements(errim) ne 0 then star.errim = imagesection(errim,p) $
     else $
        star.errim = errim(star.image, 10000., 1., 5., 0.0375, 0., 10) ;some sorta typical values...

     if total(~finite(star.image)) ne 0 then continue

     lowerlimit = (starsize-innersize)/2
     upperlimit = lowerlimit+innersize-1
     inner = (star.image)[lowerlimit:upperlimit,lowerlimit:upperlimit]
     innererr = (star.errim)[lowerlimit:upperlimit,lowerlimit:upperlimit]
     xs=findgen(innersize)+lowerlimit+0.5
     ys=findgen(innersize)+lowerlimit+0.5
     star.model = mpfit2dpeak(inner, a, xs, ys, error=innererr, /quiet)
     star.modelparams = a
     if n_elements(stars) eq 0 then stars=star else stars=[stars,star]
  endfor
  if n_elements(stars) eq 0 then return, -1
  return, stars
end
