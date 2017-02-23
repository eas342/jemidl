Function skymaskobj, obj, iband
  nxpfile = (obj->extract('filehead'))[iband]+'_NXP.fits.gz'
  
  nxp = mrdfits(nxpfile,/silent)
  mask = nxp ge 3

  mask = mask or shift( mask, 1 )
  mask = mask or shift( mask, -1 )
  mask = mask or shift( mask, 0, 1 )
  mask = mask or shift( mask, 0, -1 ) ;;extend 3x3 box
  sz = size( mask, /dim )
  
  sexcat=(obj->extract('sexcat'))[*,0]
  glpcat=obj->extract('galapagoscat')
  

  objmask0=bytarr( sz[0], sz[1] )
  for iobj=0, n_elements(sexcat)-1 do begin
     counter, iobj, n_elements(sexcat)
     glp=glpcat[iobj]
     sex=sexcat[iobj]
     xmin = glp.xmin - 50
     xmax = glp.xmax + 50
     ymin = glp.ymin - 50
     ymax = glp.ymax + 50
     
     starmult=1.
     if sex.class_star gt 0.2 then begin
        xmin -= 300
        ymin -= 300
        xmax += 300
        ymax += 300
        starmult=2.5
     endif

     xmin = xmin > 0
     ymin = ymin > 0
     xmax = xmax < (sz[0]-1)
     ymax = ymax < (sz[1]-1)
     
     dist_ellipse, dist, [xmax-xmin+1, ymax-ymin+1], $
                   sex.xwin_image-xmin, sex.ywin_image-ymin, $
                   sex.a_image/sex.b_image, sex.theta_image-90
     m1 = dist lt (15 > (sex.a_image*sex.kron_radius*1.4+30)*starmult )
     objmask0[xmin:xmax, ymin:ymax] = objmask0[xmin:xmax, ymin:ymax] or m1
  endfor

  objmask = objmask0*0B
  for yshift=-10, 10 do begin
     for xshift=-10, 10 do begin
        counter, (xshift+10+1)+21*(yshift+10), 441
        if (xshift^2 + yshift^2) gt 100 then continue
        objmask = objmask or shift( objmask0, xshift, yshift )
     endfor
  endfor

  return, byte(mask and (1-objmask))
end
