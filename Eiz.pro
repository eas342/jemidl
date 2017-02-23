Pro Eiz
  zs=findgen(100)/99.*(1.5-0.85)+0.85
  spec = bc03spec( 2.5e9, 0.017, /erg, origmass=1.e11 )
  wave = spec->wavelength()
  lumdens = spec->flux()
  dcolor=fltarr(n_elements(zs))
  for iz=0, n_elements(zs)-1 do begin
     iz0 = izmags4( wave, lumdens, zs[iz], host_ebv=0. )
     iz1 = izmags4( wave, lumdens, zs[iz], host_ebv=1./3.1 )
     color0 = iz0[0]-iz0[1]
     color1 = iz1[0]-iz1[1]
     dcolor[iz] = color1-color0
  endfor
  junk=fsc_color( /all, color=ctable )
  plot, zs, dcolor, ytitle='E(i-z)', xtitle='redshift'
  oplot, zs, dcolor/2., color=ctable.red
  oplot, zs, dcolor/4., color=ctable.blue

  legend, ['A_V=1.0', 'A_V=0.5', 'A_V=0.25'], color=[ctable.white, ctable.red, ctable.blue], $
          textcolor=[ctable.white, ctable.red, ctable.blue], charsize=2.
          
end
