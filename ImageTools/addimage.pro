;+
; NAME:
;	ADDIMAGE
;
; PURPOSE:
;       Adds one image to another with arbitrary overlap.  The size of
;       the output image is the same as the size of the first operand image.
; CALLING SEQUENCE:
;	Result = ADDIMAGE(image1, image2,...)
; INPUTS:
;       im1: Image number 1
;       im2: Image to be added to image number 1
;       w1: A coordinate in image 1 to align with w2
;       w2: A coordinate in image 2 to align with w1
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       Result: The result of image addition.  The size of output is
;       the same as the size of im1.
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       JEM, Jan, 2009. Written
;-

Function AddImage, im1, im2, w1, w2, max=max
  s1 = size(im1,/dim)
  s2 = size(im2,/dim)
  p1 = dblarr(2,2) ;;subimage coords for im1
  p2 = dblarr(2,2) ;;subimage coords for im2
  ;; ic=0 => x-coord
  ;; ic=1 => y-coord
  for ic = 0,1 do begin
     p1[ic,*] = [w1[ic]-w2[ic], w1[ic]-w2[ic]+s2[ic]-1]
     p2[ic,*] = [0,s2[ic]-1]
     if p1[ic,0] lt 0 then begin
        p2[ic,0] -= p1[ic,0]
        p1[ic,0]=0
     endif
     if p1[ic,1] ge s1[ic] then begin
        p2[ic,1] -= (p1[ic,1]-(s1[ic]-1))
        p1[ic,1] = s1[ic]-1
     endif
  endfor
  tmp = im1
  if p1[0,0] gt p1[0,1] $
     or p1[1,0] gt p1[1,1] $
     or p2[0,0] gt p2[0,1] $
     or p2[1,0] gt p2[1,1] then return, tmp
  if keyword_set(max) then $
     tmp[p1[0,0]:p1[0,1],p1[1,0]:p1[1,1]] >= im2[p2[0,0]:p2[0,1],p2[1,0]:p2[1,1]] $
  else $
     tmp[p1[0,0]:p1[0,1],p1[1,0]:p1[1,1]] += im2[p2[0,0]:p2[0,1],p2[1,0]:p2[1,1]]
  return, tmp
end
