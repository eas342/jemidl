Pro plotcore, obj, ps=ps
  resolve_obj, obj
  obj->LoadImage, band='z'
  s=obj->summary()
  ra_cluster = obj->extract('ra')
  dec_cluster = obj->extract('dec')
  id = obj->extract('clusterid')
  cd1 = linfit(s.alphawin_j2000, s.xwin_image)
  cd2 = linfit(s.deltawin_j2000, s.ywin_image)
  x_cluster = cd1[0]+cd1[1]*ra_cluster
  y_cluster = cd2[0]+cd2[1]*dec_cluster
  xmin = x_cluster-400
  xmax = x_cluster+399
  ymin = y_cluster-400
  ymax = y_cluster+399
  zim = obj->ImageSection([xmin,ymin,xmax,ymax], band='z')

  if ~keyword_set(ps) then begin
     window, 0, xsize=800, ysize=800
     tv, bytscl(zim, min=-0.01, max=0.05, top=250)
     red=fsc_color('red', 251)
     blue=fsc_color('red', 252)
     w=where( s.xmin le xmax and s.xmax ge xmin and $
              s.ymin le ymax and s.ymax ge ymin )
     for i=0, n_elements(w)-1 do begin
        npoints = 120
        phi = 2.*!pi*(findgen(npoints)/(npoints-1))
        majaxis = s[w[i]].a_image*s[w[i]].kron_radius
        minaxis = s[w[i]].b_image*s[w[i]].kron_radius
        ang = s[w[i]].theta_image*!pi/180.
        x = majaxis*cos(phi)
        y = minaxis*sin(phi)
        xprime = s[w[i]].xwin_image-xmin + x*cos(ang) - y*sin(ang)
        yprime = s[w[i]].ywin_image-ymin + x*sin(ang) + y*cos(ang)
        pts = fltarr(2, npoints)
        pts[0,*] = xprime
        pts[1,*] = yprime
        if s[w[i]].cold then color=blue else color=red
        plots, pts, /device, color=color
     endfor
  endif else begin
     set_plot, 'ps'
     device, /encapsulated, /color, bits_per_pixel=8, landscape=0
     device, xsize=8., ysize=8., /inches, xoffset=0., yoffset=0.
     device, filename=id+'core.eps'

     tvimage, bytscl(zim, min=-0.01, max=0.05, top=250)
     red=fsc_color('red', 251)
     blue=fsc_color('blue', 252)
     w=where( s.xmin le xmax and s.xmax ge xmin and $
              s.ymin le ymax and s.ymax ge ymin )
     for i=0, n_elements(w)-1 do begin
        npoints = 120
        phi = 2.*!pi*(findgen(npoints)/(npoints-1))
        majaxis = s[w[i]].a_image*s[w[i]].kron_radius
        minaxis = s[w[i]].b_image*s[w[i]].kron_radius
        ang = s[w[i]].theta_image*!pi/180.
        x = majaxis*cos(phi)
        y = minaxis*sin(phi)
        xprime = s[w[i]].xwin_image-xmin + x*cos(ang) - y*sin(ang)
        yprime = s[w[i]].ywin_image-ymin + x*sin(ang) + y*cos(ang)
        xprime /= 800
        yprime /= 800
        pts = fltarr(2, n_elements(xprime))
        pts[0,*] = xprime
        pts[1,*] = yprime
        if s[w[i]].cold then color=blue else color=red
        plots, pts, /normal, color=color
     endfor
     device, /close_file
     set_plot, 'X'
  endelse

  obj->FreeImages
  
end
