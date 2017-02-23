;+
; NAME:
;	PSFRESID
;
; PURPOSE:
;       Corrects an oversampled psf by adding star fit residuals to it.
; CALLING SEQUENCE:
;	psf = PSFRESID( stars, oversample_factor, overpsf)
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
;       overpsf: the oversampled psf to fit to the stars.
;       oversample_factor: Factor by which the PSF has been
;                          oversampled
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       Result:  The median residual image after fitting the PSF to stars.
;
; EXAMPLE:
;
; NOTES:
;       Josh coords place the origin at the lower left corner of the
;       lower left pixel.
;
; MODIFICATION HISTORY:
;       JEM, Jan, 2009. Written
;-

Function PSFresid, stars, oversample_factor, overpsf, _extra=extra
; extra: /interp
  oversize = (size(overpsf,/dim))[0]
  if size(stars,/tname) ne 'STRUCT' then return, 0
  coverresid = dblarr(oversize,oversize,n_elements(stars))*!values.f_nan ;; residuals for each star
  for istar=0, n_elements(stars)-1 do begin
     a=fitstar(stars[istar],overPSF,oversample_factor, _extra=extra)
     residsize = (size(a.resid,/dim))[0]*oversample_factor
     overresid1 = oversampleimage( a.resid, oversample_factor, _extra=extra )/a.acoeff[0]  ;;acoeff needed to normalize...
     xlow=oversample_factor-a.xshift
     ylow=oversample_factor-a.yshift
     if xlow lt 0 $
        or ylow lt 0 $
        or xlow gt 2*oversample_factor $
        or ylow gt 2*oversample_factor then begin ;;ignore stars where the fit is very far from center of postage stamp
        coverresid[*,*,istar]=!values.f_nan
        continue
     endif
     coverresid[xlow:xlow+residsize-1, ylow:ylow+residsize-1, istar]=overresid1
  endfor
  return, median( coverresid, dim=3 ) ;;smash using median
end
