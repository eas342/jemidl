Pro zhist, objs

  colors = Obj_New("IDLgrPalette")
  colors->LoadCT, 3
  colors->GetProperty, Red=r, Green=g, Blue=b
  Obj_Destroy, colors
  TVLCT, r, g, b
  red = getcolor('red', !d.table_size-2)
  green = getcolor('green', !d.table_size-3)

  thisDevice=!d.name
  set_plot, 'ps'
  device, /encapsulated, color=1, bits_per_pixel=8
  device, xsize=8., ysize=10., /inches
  device, filename='zhist.eps'

  !p.multi=[0,1,25,0,1]
  for iobj=0, n_elements(objs)-1 do begin
     position = [0.1, 1.0-(iobj+2.d)/(n_elements(objs)+3), 0.9, 1.0-(iobj+1.d)/(n_elements(objs)+3)]
     speccat = objs[iobj]->speccat()
     if size(speccat,/tname) ne 'STRUCT' then continue
     w=where(speccat.zqual ne 'N' and speccat.zqual ne 'F')
     if w[0] eq -1 then continue
     objs[iobj]->GetProperty, clustername=clustername
     case clustername of 
        'A' : w1=where(speccat.z ge 1.44 and speccat.z le 1.46)
        'B' : w1=where(speccat.z ge 1.10 and speccat.z le 1.14)
        'C' : w1=where(speccat.z ge 0.94 and speccat.z le 1.00)
        'E' : w1=where(speccat.z ge 1.00 and speccat.z le 1.05)
        'F' : w1=where(speccat.z ge 1.10 and speccat.z le 1.15)
        'H' : w1=where(speccat.z ge 1.225 and speccat.z le 1.235)
        'K' : w1=where(speccat.z ge 1.40 and speccat.z le 1.42)
        'N' : w1=where(speccat.z ge 1.01 and speccat.z le 1.03)
        'O' : w1=where(speccat.z ge 1.017 and speccat.z le 1.03)
        'R' : w1=where(speccat.z ge 1.20 and speccat.z le 1.22)
        'T' : w1=where(speccat.z ge 0.97 and speccat.z le 0.98)
        'U' : w1=where(speccat.z ge 1.02 and speccat.z le 1.045)
        'V' : w1=where(speccat.z ge 0.90 and speccat.z le 0.91)
        'X' : w1=where(speccat.z ge 1.05 and speccat.z le 1.15)
        'Y' : w1=where(speccat.z ge 1.22 and speccat.z le 1.25)
        'Z' : w1=where(speccat.z ge 1.39 and speccat.z le 1.40)
        else : w1=0
     endcase
     if n_elements(w1) ge 5 then $
        plotzhist, speccat.z, w1, bin=0.005, xrange=[0.85, 1.5], yrange=[0,15], pos=position, xstyle=1, ystyle=1, gcolor=red $
     else $
        plothist, speccat.z, bin=0.005, xrange=[0.85, 1.5], yrange=[0,15], pos=position, xstyle=1, ystyle=1
     print, clustername
     xyouts, position[0]-0.03, position[1]+0.01, clustername, /normal
  endfor
  device, /close_file
  set_plot, thisDevice  
end
