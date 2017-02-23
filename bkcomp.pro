Pro bkcomp, obj, bk, bmagauto=bmagauto, jmagauto=jmagauto, $
            bcolor=bcolor, jcolor=jcolor, $
            bre=bre, jre=jre, $
            bgini=bgini, jgini=jgini, $
            basym=basym, jasym=jasym, $
            ps=ps
  nbk=n_elements(bk)
  s=obj->summary()
  c=obj->color(/clean)
  ebv=obj->extract('ebv')
  A_i = 1.973*ebv
  A_z = 1.472*ebv

  s.zmag_auto -= A_z
  c -= (A_i-A_z)
  s.zmag_auto -= (s.zmag_auto-21)*0.1/3.+0.2

  bmagauto=fltarr(nbk)
  jmagauto=fltarr(nbk)
  bcolor=fltarr(nbk)
  jcolor=fltarr(nbk)
  bre=fltarr(nbk)
  jre=fltarr(nbk)
  bgini=fltarr(nbk)
  jgini=fltarr(nbk)
  basym=fltarr(nbk)
  jasym=fltarr(nbk)
  for ibk=0, nbk-1 do begin
     b=bk[ibk]
     w=obj->nearestobjs(b.alpha_j2000, b.delta_j2000)

     bmagauto[ibk]=b.mag_auto
     jmagauto[ibk]=s[w[0]].zmag_auto

     bcolor[ibk]=b.iz
     jcolor[ibk]=c[w[0]]

     bre[ibk]=b.re*0.03
     jre[ibk]=s[w[0]].re*0.05

     bgini[ibk]=b.gini
     jgini[ibk]=s[w[0]].gini2

     basym[ibk]=b.asym
     jasym[ibk]=s[w[0]].asym2

  endfor
  if keyword_set(ps) then begin
     thisDevice=!d.name
     set_plot, 'ps'
     device, /encapsulated, /color, bits_per_pixel=8, landscape=0
     device, xsize=8., ysize=10., /inches
     device, filename = obj->extract('clustername')+'.eps'
     charsize=1.5
  endif else begin
     window,0, xsize=1000, ysize=1000
     charsize=2.0
  endelse

  !p.multi=[0,2,5]
  !y.omargin = [2,6]

  plot, jmagauto, bmagauto, ps=3, xrange=[19,25], yrange=[19,25], title='MAG_AUTO', charsize=charsize, $
        xtitle='josh', ytitle='ben'
  oplot, [19,25], [19,25]
  plot, jmagauto, bmagauto-jmagauto, ps=3, xrange=[19,25], yrange=[-1,1], $
        title=textoidl('\DeltaMAGAUTO'), charsize=charsize
  oplot, [19,25], [0,0]
  junk = biweight_mean(bmagauto-jmagauto,magsig)
  print, 'magsig: ', magsig
  xyouts, 19, 0.5, textoidl('\sigma=')+string(magsig, format='(F5.3)')

  plot, jcolor, bcolor, ps=3, xrange=[-0.5,1.5], yrange=[-0.5,1.5], title='color', charsize=charsize, $
        xtitle='josh', ytitle='ben'
  oplot, [-1,2], [-1,2]
  plot, jmagauto, bcolor-jcolor, ps=3, xrange=[19,25], yrange=[-0.5,0.5], $
        title=textoidl('\Deltacolor'), charsize=charsize
  oplot, [19,25], [0,0]
  wc = where(finite(bcolor-jcolor))
  junk = biweight_mean((bcolor-jcolor)[wc],colorsig)
  print, 'colorsig: ', colorsig
  xyouts, 19, 0.4, textoidl('\sigma=')+string(colorsig, format='(F5.3)')

  plot, jre, bre, ps=3, xrange=[0,1], yrange=[0,1], title='Re[arcsec]', charsize=charsize, $
        xtitle='josh', ytitle='ben'
  oplot, [0,2], [0,2]
  plot, jmagauto, bre-jre, ps=3, xrange=[19,25], yrange=[-0.5,0.5], $
        title=textoidl('\DeltaRe[arcsec]'), charsize=charsize
  oplot, [19,25], [0,0]
  junk = biweight_mean(bre-jre, resig)
  print, 'resig: ', resig
  xyouts, 19, 0.4, textoidl('\sigma=')+string(resig, format='(F5.3)')

  plot, jgini, bgini, ps=3, xrange=[0.2,0.8], yrange=[0.2,0.8], title='gini', charsize=charsize, $
        xtitle='josh', ytitle='ben'
  oplot, [0,1], [0,1]
  plot, jmagauto, bgini-jgini, ps=3, xrange=[19,25], yrange=[-0.2, 0.2], $
        title=textoidl('\Deltagini'), charsize=charsize
  oplot, [19,25], [0,0]
  w=where(finite(bgini-jgini))
  junk = biweight_mean((bgini-jgini)[w], ginisig)
  print, 'ginisig: ', ginisig
  xyouts, 19, 0.13, textoidl('\sigma=')+string(ginisig, format='(F5.3)')

  plot, jasym, basym, ps=3, xrange=[0,0.5], yrange=[0,0.5], title='asym', charsize=charsize, $
        xtitle='josh', ytitle='ben'
  oplot, [0,1], [0,1]
  plot, jmagauto, basym-jasym, ps=3, xrange=[19,25], yrange=[-0.3, 0.3], $
        title=textoidl('\Deltaasym'), charsize=charsize
  oplot, [19,25], [0,0]
  w=where(finite(basym-jasym))
  junk = biweight_mean((basym-jasym)[w], asymsig)
  print, 'asymsig: ', asymsig
  xyouts, 19, 0.25, textoidl('\sigma=')+string(asymsig, format='(F5.3)')

  xyouts, 0.5, 0.95, /normal, obj->extract('clustername'), align=0.5, charsize=charsize*1.5
  !p.multi=0

  if keyword_set(ps) then begin
     device, /close_file
     set_plot, thisDevice
  endif
  
end
