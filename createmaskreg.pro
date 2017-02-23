Pro CreateMaskReg, objfile, file1, file2, outfile, widthin=widthin
  readcol, objfile, name, rah, ram, ras, ddeg, dmin, dsec, format='A,X,X,I,I,F,I,I,F,X,X,X,X'
  readcol, file1, x_star, y_star, min_y, max_y, percent, name2, mag, format='F,F,F,F,F,X,X,X,A,F'
  readcol, file2, param, param2, default, pa, format='A,A,F,F'
  pa = (pa[where(param eq 'Position' and param2 eq 'angle:')])[0]
  
  openw, lun, outfile, /get_lun
  for iname=0, n_elements(name2)-1 do begin
     if name2[iname] eq 'CENTER' then continue
     if strpos(name2[iname], 'star') ne -1 then begin
        width=4.
        length=4.
        angle=pa
        color='red'
        w2 = where(name2[iname] eq name)
        rah1 = rah[w2]
        ram1 = ram[w2]
        ras1 = ras[w2]
        ddeg1 = ddeg[w2]
        dmin1 = dmin[w2]
        dsec1 = dsec[w2]
     endif else begin
        w2 = where(name2[iname] eq name)
        rah1 = rah[w2]
        ram1 = ram[w2]
        ras1 = ras[w2]
        ddeg1 = ddeg[w2]
        dmin1 = dmin[w2]
        dsec1 = dsec[w2]
        if n_elements(widthin) eq 0 then width=1. $
        else width=widthin
        length=1.
        angle=pa
        color='green'
        str=string(format='(%"box(%02i:%02i:%06.3f, %+03i:%02i:%06.3f, %6.3f\", %8.3f\", %+8.3f) # color=%s")', $
                   rah1, ram1, ras1, ddeg1, dmin1, dsec1, width, length, angle, color)
        printf, lun, str

        if n_elements(widthin) eq 0 then width=1. $
        else width=widthin
        length=(min_y[iname]-max_y[iname])/0.725
        angle=pa
        color='blue'
        shift = length*(percent[iname]/100. - 0.5)
        shiftalpha = shift*sin(-pa*!dpi/180.d)
        shiftdelta = -shift*cos(-pa*!dpi/180.d)

        w2 = where(name2[iname] eq name)
        rah1 = rah[w2]
        ram1 = ram[w2]
        ras1 = ras[w2]
        ddeg1 = ddeg[w2]
        dmin1 = dmin[w2]
        dsec1 = dsec[w2]
        ra = 15.*ten(rah[w2], ram[w2], ras[w2])
        dec = ten(ddeg[w2], dmin[w2], dsec[w2])
        ra += shiftalpha/cos(dec*!dpi/180.d)/3600.
        dec += shiftdelta/3600.
        ra=sixty(ra/15.)
        dec=sixty(dec)
        rah1=ra[0]
        ram1=ra[1]
        ras1=ra[2]
        ddeg1=dec[0]
        dmin1=dec[1]
        dsec1=dec[2]
     endelse
     str=string(format='(%"box(%02i:%02i:%06.3f, %+03i:%02i:%06.3f, %8.3f\", %8.3f\", %+8.3f) # text={%s} color=%s")', $
                rah1, ram1, ras1, ddeg1, dmin1, dsec1, width, length, angle, name2[iname], color)
     print, str
     printf, lun, str
  endfor
  close, lun
  free_lun, lun
end
