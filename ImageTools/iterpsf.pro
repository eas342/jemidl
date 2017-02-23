;+
; NAME:
;	ITERPSF
;
; PURPOSE:
;       Iteratively construct a psf from stars.
; CALLING SEQUENCE:
;	overpsf = ITERPSF( stars, niter, oversample_factor,
;	                   [resids=resids, overpsfs=overpsfs] )
; INPUTS:
;       stars: an array of star structures with at least the following
;              elements:
;              star.image: an NxN image of a star
;              star.modelparams: an array of gaussian fit
;                    coefficients:
;                            modelparams[0] : an offset
;                            modelparams[1] : amplitude
;                            modelparams[2] : x-sigma
;                            modelparams[3] : y-sigma
;                            modelparams[4] : x-center
;                            modelparams[5] : y-center
;       model = model[0] + model[1]*exp(-0.5((x-model[4])/model[2])^2+ ...)
;       niter: number of residual correcting iterations to perform
;       oversample_factor: factor by which psf is oversampled
;
; KEYWORD PARAMETERS:
;
; KEYWORD OUTPUT:
;       resids: the oversampled median residual after each iteration
;       overpsfs: the oversampled psf after each iteration
;
; OUTPUTS:
;       overpsf: The final oversampled PSF.
; EXAMPLE:
;
; NOTES:
;       Josh coords place the origin at the lower left corner of the
;       lower left pixel.
;
; MODIFICATION HISTORY:
;       JEM, Jan, 2009. Written
;-

Function IterPSF, stars, niter, $
                  oversample_factor, $
                  resids=resids, $
                  overPSFs=overPSFs, $
                  _extra=extra
; extras: psfsize=psfsize
;         /plot
;         /interp
  junk = getPSF( stars, oversample_factor=oversample_factor, $
                 overPSF=overPSF0, /binfrac, /interp, _extra=extra )
  if n_elements(junk) eq 1 then begin
     message, 'error with PSF'
     return, 0
  endif
  oversize=(size( overPSF0, /dim ))[0]
  overPSFs=dblarr( oversize, oversize, niter+1 )
  resids=dblarr( oversize, oversize, niter )
  overPSFs[*,*,0]=overPSF0
  for i=0, niter-1 do begin
     overPSF = overPSFs[*,*,i]
     resid = PSFresid( stars, oversample_factor, overPSF, _extra=extra )
     resids[*,*,i]=resid
     overPSFs[*,*,i+1]=overPSFs[*,*,i]+resid
  endfor
;;trim edges which may not have been evenly sampled.
  overPSFs = overPSFs[oversample_factor:oversize-oversample_factor-1, $
                      oversample_factor:oversize-oversample_factor-1, $
                      *]
  resids = resids[oversample_factor:oversize-oversample_factor-1, $
                  oversample_factor:oversize-oversample_factor-1, $
                  *]
  return, overPSFs[*,*,niter]
end
