Function Galapagos1, sexcat
  ;; setup catalog to fill in next
  galapagoscat = { galapagoscat, $
                   clusterid:  '',$
                   galid:      0L, $
                   xmin:       0L, $
                   xmax:       0L, $
                   ymin:       0L, $
                   ymax:       0L, $
                   simfit:     ptr_new(), $
                   maskfit:    ptr_new(), $
                   simfit2:    ptr_new(), $
                   maskfit2:   ptr_new(), $
                   majaxis:    0.d }

  nobjs = (size(sexcat,/dim))[0]
  galapagoscat = replicate( galapagoscat, nobjs )
  sex0 = sexcat[*,0]

  theta = sex0.theta_image*!dpi/180.d
  sth = abs( sin(theta) )
  cth = abs( cos(theta) )
  
  ;; size of postage stamp image
  xsize = 3.0*sex0.a_image*sex0.kron_radius*(cth + (1. - sex0.ellipticity)*sth)
  ysize = 3.0*sex0.a_image*sex0.kron_radius*(sth + (1. - sex0.ellipticity)*cth)
  
  xmin = long( sex0.xwin_image-1 - xsize/2. )
  xmax = long( sex0.xwin_image-1 + xsize/2. ) + 1
  ymin = long( sex0.ywin_image-1 - ysize/2. )
  ymax = long( sex0.ywin_image-1 + ysize/2. ) + 1
  
  ;; size of ellipses for each object
  majaxis = sex0.a_image*sex0.kron_radius*1.5

  ;; loop through each object in sexcat
  for i=0, n_elements(galapagoscat)-1 do begin
     counter, i+1, n_elements(galapagoscat), 'Executing GALAPAGOS on object '
     sex1=sex0[i]
     nx = xmax[i] - xmin[i] + 1
     ny = ymax[i] - ymin[i] + 1
     dist_ellipse, imellipse, [nx, ny], $
                   sex1.xwin_image - 1 - xmin[i], sex1.ywin_image - 1 - ymin[i], $
                   sex1.a_image/sex1.b_image, $
                   sex1.theta_image - 90.
     orig_ellipse = imellipse lt majaxis[i]
     orig_ellipse2 = imellipse lt majaxis[i]/1.5  ;; for less conservative 2nd attempt GALFITing
     
     ;; reset simfit and maskfit for each galaxy
     delvarx, simfit, maskfit, simfit2, maskfit2
     ;; find other objects whose postage stamps overlap with
     ;; this object's postage stamp.
     w = where( xmin le xmax[i] and xmax ge xmin[i] and $
                ymin le ymax[i] and ymax ge ymin[i] )
     if w[0] ne -1 then begin
        for k=0, n_elements(w)-1 do begin
           j = w[k]
           maskgal=sex0[j]
           if j eq i then continue  ;; don't mask target object!
           dist_ellipse, ellipse, [nx, ny], $
                         maskgal.xwin_image - xmin[i] - 1, maskgal.ywin_image - ymin[i] - 1, $
                         maskgal.a_image/maskgal.b_image, $
                         maskgal.theta_image - 90.
           e = ellipse lt majaxis[j]
           e2 = ellipse lt (majaxis[j]/1.5)
           if total(e and orig_ellipse) gt 0 then begin  ;; if ellipses intersect
              if n_elements(simfit) eq 0 then simfit = j $
              else simfit = [simfit, j]
           endif else if total(e) gt 0 then begin  ;; if mask ellipse intersects postage stamp
              if n_elements(maskfit) eq 0 then maskfit = j $
              else maskfit = [maskfit, j]
           endif
           ;; fill in less conservative mask/sim fit
           if total(e2 and orig_ellipse2) gt 0 then begin  ;; ellipses intersect
              if n_elements(simfit2) eq 0 then simfit2 = j $
              else simfit2 = [simfit2, j]
           endif else if total(e2) gt 0 then begin  ;; postage stamp intersects
              if n_elements(maskfit2) eq 0 then maskfit2 = j $
              else maskfit2 = [maskfit2, j]
           endif
        endfor
     endif
     if n_elements(simfit) ne 0 then galapagoscat[i].simfit = Ptr_New(simfit)
     if n_elements(maskfit) ne 0 then galapagoscat[i].maskfit = Ptr_New(maskfit)
     if n_elements(simfit2) ne 0 then galapagoscat[i].simfit2 = Ptr_New(simfit2)
     if n_elements(maskfit2) ne 0 then galapagoscat[i].maskfit2 = Ptr_New(maskfit2)
  endfor
  print

  galapagoscat.clusterid = sex0.clusterid
  galapagoscat.galid = sex0.galid
  galapagoscat.xmin = xmin
  galapagoscat.xmax = xmax
  galapagoscat.ymin = ymin
  galapagoscat.ymax = ymax
  galapagoscat.majaxis = majaxis
  return, galapagoscat
end
