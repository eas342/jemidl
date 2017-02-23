Function likebkg, data, bkgparams
  return,  bkgparams.d $
           + bkgparams.e*data.mag $
           + bkgparams.f*data.color $
           + bkgparams.g*data.mag^2
end

Function pq_limits, x
  common JEM$_bkg, magrange, colorrange, bkgparams, data, gals ;;colorrange
  return, [colorrange[0], colorrange[1]]
end

;;scalar mag, vector color
Function fieldintegrand, mag, color
  common JEM$_bkg ;; bkgparams
  return, 159.28*likebkg( {mag:replicate(mag, n_elements(color)), $
                           color:color}, $
                          bkgparams )
end

Function loglikebkg, abscissa
  common JEM$_bkg ;;data, magrange
  d=abscissa[0]
  e=abscissa[1]
  f=abscissa[2]
  g=abscissa[3]
  bkgparams = { d:d, e:e, f:f, g:g }
  integral = int_2d( 'fieldintegrand', magrange, 'pq_limits', 20 )
  out = total(alog(likebkg(data, bkgparams)))-integral
 return, -out
end


Pro bkg, start, scale, magrange=magrange1, colorrange=colorrange1
  common JEM$_bkg
  if n_elements(gals) eq 0 then begin
     restore, '/home/scpdata03/goods/cmd/cmds.sav'
     gals = *(*cmds[0]).gals
     for i=1, n_elements(cmds)-1 do begin
        gals=[gals, *(*cmds[i]).gals]
     endfor
  endif
  data = {mag:0., color:0.}
  data = replicate( data, n_elements(gals) )
  data.mag = gals.zmag
  data.color = gals.iz
  if n_elements(magrange1) eq 0 then magrange = [20., 24.] $
  else magrange=magrange1
  if n_elements(colorrange1) eq 0 then colorrange = [0.6, 1.2] $
  else colorrange=colorrange1
  w=where( data.mag gt magrange[0] $
           and data.mag le magrange[1] $
           and data.color gt colorrange[0] $
           and data.color le colorrange[1] )
  data = data[w]

  if n_elements(start) eq 0 then start = [-1923.81, 170.483, -762.931, -1.82123]/159.
  if n_elements(scale) eq 0 then scale = [10., 1., 10., 0.01]

  result = downhillsimplex( 'loglikebkg', start, scale, $
                            tol=0.00001, fval=fval, miniter=500, maxiter=20000 )
  print
  print, result
  print
  print, fval[0]
  print, n_elements(data)
  bkgparams = { d:result[0], e:result[1], f:result[2], g:result[3] }
  print, int_2d( 'fieldintegrand', magrange, 'pq_limits', 20 )
end
