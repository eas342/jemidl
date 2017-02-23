;+
; NAME:
;	LOWPASS
;
; PURPOSE:
;       This function returns a filtered waveform computed by zeroing
;       the high frequency components of a discrete Fourier transform.
; CALLING SEQUENCE:
;	Result = LOWPASS(data, ...)
; INPUTS:
;       data: A vector of waveform values.
;       threshold: A scalar below which frequency components are preserved.
; KEYWORD PARAMETERS:
;       f: An optional output of the frequency values derived from the
;       input vector.
; OUTPUTS:
;       Result: The filtered waveform vector.
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       JEM, July, 2008. Written
;-

Function Lowpass, data, threshold, f=f
  a=fft(data)
  n=n_elements(data)
  n21 = n/2+1
  f = dindgen(n)
  if n_elements(data) mod 2 eq 1 then $
     f[n21] = n21 - n + findgen(n21-1) $
  else $
     f[n21] = n21 - n + findgen(n21-2)
  f=f/n
  w=where(abs(f) GT threshold)
  if w[0] ne -1 then a[w]=0.
  return, real_part(fft(a, /inverse))
end
