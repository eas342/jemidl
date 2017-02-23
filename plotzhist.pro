;+
; NAME:
;	PLOTZHIST
;
; PURPOSE:
;       Plot a histogram of redshifts and optionally fit
;       and overplot a gaussian derived from the biweight estimator.
; CALLING SEQUENCE:
;	PLOTZHIST, z, w, bin=bin, zmean=zmean, sigma=sigma, gcolor=gcolor
; INPUTS:
;       z: array of redshifts.
;       w: indices of z with which to fit a guassian
;
; KEYWORD PARAMETERS:
;       bin: histogram binning resolution
;       gcolor: color to with which to overplot gaussian
;       fcolor: fill color
;
; OUTPUTS:
;       zmean: result of biweight fit
;       sigma: the fit dispersion in units of redshift
;       
; EXAMPLE:
;
; z = [0.1, 0.09, 0.09, 0.092, 0.094, 0.11, 0.1, 0.08, 0.01, 0.02]
; w = where(z gt 0.05)
; plotzhist, z, w, zmean=zmean, sigma=sigma 
; print, zmean, sigma*299792/(1.+zmean)
;
; MODIFICATION HISTORY:
;       JEM, July, 2008. Written
;-

Pro PlotZHist, z, w, bin=bin, zmean=zmean, sigma=sigma, fcolor=fcolor, gcolor=gcolor, _extra=extra
  ;; defaults...
  if n_elements(bin) eq 0 then bin=0.05
  if n_elements(gcolor) eq 0 then gcolor = !d.table_size/2
  plothist, z, bin=bin, _extra=extra
  if n_elements(w) lt 2 then return
  plothist, z[w], bin=bin, _extra=extra, /fill, /over, fcolor=fcolor
  if n_elements(w) lt 5 then return
  ;; do the fit
  zmean=biweight_mean(z[w], sigma)
  ;; overplot up to z=2.0
  x=findgen(5000)/4999.*2
  ;; make gaussian have the same area as z[w] in the histogram
  y=exp(-(x-zmean)^2/(2*sigma^2))*n_elements(w)*bin/sqrt(2.*!dpi)/sigma
  w1=where(y ge max(y)*0.001)
  oplot, x[w1], y[w1], color=gcolor, _extra=extra
end
