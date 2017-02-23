Function objmask, sexcat, xsize, ysize
  mask = bytarr( xsize, ysize )
  for iobj=0, n_elements(sexcat)-1 do begin
     counter, iobj+1, n_elements(sexcat)
     sex = sexcat[iobj]
     maxsize = sex.a_image*sex.kron_radius*2.0
     xmin = long(sex.xwin_image - maxsize-1)
     xmax = long(sex.xwin_image + maxsize)
     ymin = long(sex.ywin_image - maxsize-1)
     ymax = long(sex.ywin_image + maxsize)
     starmult=1.
     if sex.class_star gt 0.2 then begin ;;stars get extra masking...
        xmin -= 300
        xmax += 300
        ymin -= 300
        ymax += 300
        starmult=2.0
     endif
     xmin = xmin > 0
     ymin = ymin > 0
     xmax = xmax < (xsize-1)
     ymax = ymax < (ysize-1)

     if xmax lt xmin then continue
     if ymax lt ymin then continue

     dist_ellipse, dist, [xmax-xmin+1, ymax-ymin+1], $
                   sex.xwin_image-xmin, sex.ywin_image-ymin, $
                   sex.a_image/sex.b_image, sex.theta_image-90.
     m1 = dist lt (10 > (maxsize))*starmult
     mask[xmin:xmax, ymin:ymax] = mask[xmin:xmax, ymin:ymax] or m1
  endfor
  return, mask
end
