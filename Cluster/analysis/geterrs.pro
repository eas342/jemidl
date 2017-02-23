;repeated measurements errorbars...
Function geterrs, obj, _extra=extra
  readcol, '/home/scpdata02/photerror/errmatrix.dat', format='A,F,F,F,F,F,F,F,F', $
           name, col1, col2, col3, col4, col5, col6, col7, col8, /silent
  data = transpose([[col1], [col2], [col3], [col4], [col5], [col6], [col7], [col8]])
  mags = data[*,0]
  data = data[*,indgen(50)*2+1]
  name = name[indgen(50)*2+1]
  s=obj->summary()
  c=obj->color(_extra=extra, /silent)
  m=s.zmag_auto
  w=where(strmid(name,0,1) eq obj->extract('clusterid'))
  ierr = interpol(data[*,w[0]], mags, m)
  zerr = interpol(data[*,w[1]], mags, m+c)
  return, sqrt(ierr^2+zerr^2)
end
