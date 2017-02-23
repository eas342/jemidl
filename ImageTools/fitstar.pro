;+
; NAME:
;	FITSTAR
;
; PURPOSE:
;       Fit an oversampled PSF to a star
; CALLING SEQUENCE:
;	result = FITSTAR( star, overPSF, oversample_factor,
;	[plot=plot] )
;
; INPUTS:
;       star: a structure with at least the following fields:
;              star.image: an NxN image of a star
;              star.errim: an NxN error image of the star
;              star.modelparams: an array of gaussian fit
;                    coefficients:
;                            modelparams[0] : an offset
;                            modelparams[1] : amplitude
;                            modelparams[2] : x-sigma
;                            modelparams[3] : y-sigma
;                            modelparams[4] : x-center
;                            modelparams[5] : y-center
;       overPSF: oversampled PSF
;       oversample_factor: factor by which psf is oversampled
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;       result: The fitting results for the star.  Fields are:
;            image: image of the star
;            model: fit model
;            resid: residual
;            fracresid55: the mean fractional residual in the central
;                         5x5 pixel region
;            totals: the flux in the image and the model in the
;                    central 5x5 region
;            xshift: the number of pixels shifted in the x-direction
;            yshift: the number of pixels shifted in the y-direction
;            acoeff: the best fit PSF amplitude
;            chi2perdof: the reduced chi-squared of the fit.

; EXAMPLE:
;
; NOTES:
;       Josh coords place the origin at the lower left corner of the
;       lower left pixel.
;
; MODIFICATION HISTORY:
;       JEM, Jan, 2009. Written
;-

Function modelstar, X, P, $
                    overPSF=overPSF, $
                    oversample_factor=oversample_factor, $
                    acoeff=acoeff, $
                    data=data, $
                    ierr=ierr, $
                    chisq=chisq
  imagesize = sqrt(n_elements(data))
  oversize = (size(overPSF,/dim))[0]

  i   = round(p[0]) ;;x-shift
  j   = round(p[1]) ;;y-shift
  tmpPSF = shift(overPSF, i, j)
  overlow = oversample_factor
  overhigh = oversize-oversample_factor-1
;; full-pixel scale shifted PSF estimate
  model = reform(undersampleimage(tmpPSF[overlow:overhigh, overlow:overhigh],oversample_factor), imagesize*imagesize)
  const = dblarr(imagesize*imagesize)+1.d  ;; fit a linear combination of constant and PSF estimate
  chisq = computechi2(data, ierr, [[model], [const]], acoeff=acoeff, yfit=yfit)
  model = reform(yfit,imagesize,imagesize) ;; best-fit linear combination
  return, model
end

Function fitstar, star, overPSF, oversample_factor, plot=plot
  if n_elements(oversample_factor) eq 0 then oversample_factor=9
  starsize = (size( star.image, /dim ))[0]
  oversize = (size( overPSF, /dim ))[0]
  imagesize = oversize/oversample_factor-2 ;;the -2 is because we need room to shift the psf around...

  starlow = (starsize-imagesize)/2
  starhigh = starlow+imagesize-1

  image = star.image[starlow:starhigh,starlow:starhigh]
  data = reform( image, imagesize*imagesize )
  err = star.errim[starlow:starhigh,starlow:starhigh]
  ierr  = reform( 1.d/err, imagesize*imagesize )
  w=where(1-finite(ierr))
  if w[0] ne -1 then ierr[w]=0.d  ;;force NaNs -> 0.0
  x = indgen(n_elements(image)) ;; dummy variable to make mpfitfun happy

  p0=[1.d,1.d]
  parinfo1 = { value:0.d, $
               step:4.d, $
               mpminstep:0.5d, $
               mpmaxstep:10.d, $
               parname:'xshift' }
  parinfo2 = { value:0.d, $
               step:4.d, $
               mpminstep:0.5d, $
               mpmaxstep:10.d, $
               parname:'yshift' }
  parinfo = [parinfo1, parinfo2]
  results = mpfitfun( 'modelstar', x, image, err, $
                      functargs={overPSF:overPSF, data:data, ierr:ierr, oversample_factor:oversample_factor}, $
                      parinfo=parinfo, dof=dof, $
                      /quiet )
  results=round(results)  ;; best-fit x-shift and y-shift for oversampled PSF estimate

  ;;now check surrounding 8 subpixels...  I forgot why I do this.
  chi2=dblarr(3,3)
  for i=-1,1 do begin
     for j=-1,1 do begin
        model = modelstar( x, results+[i,j], overPSF=overPSF, data=data, ierr=ierr, chisq=chisq, oversample_factor=oversample_factor )
        chi2[i+1,j+1]=chisq
     endfor
  endfor
  min=min( chi2, m )  ;; pick the best one.
  ind=array_indices( chi2, m )
  results += [-1,-1]+ind


  model = modelstar( x, results, overPSF=overPSF, data=data, ierr=ierr, acoeff=acoeff, oversample_factor=oversample_factor, chisq=chisq )
  resid = image-model ;; residual of best-fit shifted oversampled PSF + constant model
  if keyword_set(plot) then begin
     atv, [image,model,resid]
  endif

  in5low = (imagesize-5)/2 ;;forgot why I collected these too...
  in5high = in5low+5-1
  out = { image:image, $
          model:model, $
          resid:resid, $
          fracresid55: mean(abs(((image-model)/image)[in5low:in5high,in5low:in5high])), $
          totals:[total(image[in5low:in5high,in5low:in5high]), $
                  total(model[in5low:in5high,in5low:in5high])], $
          xshift: results[0], $
          yshift: results[1], $
          acoeff: acoeff, $
          chi2perdof: chisq/dof }

  return, out
end
