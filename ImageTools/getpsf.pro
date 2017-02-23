;+
; NAME:
;	GETPSF
;
; PURPOSE:
;       Generate a PSF from images of stars.
; CALLING SEQUENCE:
;	Result = GETPSF(stars, ...)
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
;
; KEYWORD PARAMETERS:
;       binfrac: use only the central subpixels for each oversampled
;                pixel
;       interp:  use the cubic interpolation feature of the IDL
;                function congrid with oversampling.
; OUTPUTS:
;       Result:  The 1-sampled PSF.
;       overpsf: The oversampled PSF.
; EXAMPLE:
;
; NOTES:
;       Josh coords place the origin at the lower left corner of the
;       lower left pixel.
;
; MODIFICATION HISTORY:
;       JEM, Jan, 2009. Written
;-

Function GetPSF, stars, $
                 psfsize=psfsize, $
                 oversample_factor=oversample_factor, $
                 overpsf=overpsf, $
                 _extra=extra
;extras...       binfrac=binfrac,
;                interp=interp

  if n_elements(psfsize) eq 0 then psfsize=33
  if n_elements(oversample_factor) eq 0 then oversample_factor = 9
  oversize = psfsize*oversample_factor
  if size(stars, /tname) ne 'STRUCT' then begin
     message, 'error with stars structure'
     return, 0
  endif

  starsize=(size(stars[0].image,/dimen))[0]
  cpsf = dblarr( oversize, oversize, n_elements(stars) )
  cpsf *= !values.f_nan

  halfoversize = (oversize-1)/2
  lowlimit = (oversize-1)/2
  highlimit = starsize*oversample_factor-1-lowlimit

  for istar=0, n_elements(stars)-1 do begin
     xc = fix(stars[istar].modelparams[4]*oversample_factor)
     yc = fix(stars[istar].modelparams[5]*oversample_factor)
     if (xc lt lowlimit or xc gt highlimit or $
         yc lt lowlimit or yc gt highlimit) then continue
     tmp = oversampleimage( stars[istar].image, oversample_factor, _extra=extra )
     cpsf[*,*,istar] = tmp[xc-halfoversize:xc+halfoversize, $  ;;extract centered subpixel array
                           yc-halfoversize:yc+halfoversize]
     low = halfoversize-oversample_factor*3
     high= halfoversize+oversample_factor*3-1
     cpsf[*,*,istar] /= total(cpsf[low:high,low:high,istar], /nan)  ;;normalize by central 3x3 full-pixel flux
  endfor
  overpsf = dblarr( oversize, oversize )
  overpsf = median( cpsf, dim=3 )  ;;oversampled PSF estimate
  cpsfsmash = undersampleimage( overpsf, oversample_factor )  ;;full-pixel sampled PSF estimate
  overpsf /= max(overpsf)  ;;normalize
  return, cpsfsmash/max(cpsfsmash)
end
