;repeated measurements errorbars...
Function geterrs2, obj, _extra=extra
  id=obj->extract('clusterid')
  restore, '/home/scpdata02/photerrorxpsf/'+id+'/'+id+'errtable.sav'
;  x=errtable.zmag_auto
  s=obj->summary()
  c=obj->color(/xpsf)
  x=s.zmag_auto
  zerr = p[0]+p[1]*(x-20)+p[2]*(x-20)^2 > errmin
  w=where(x le 20)
  if w[0] ne -1 then zerr[w] = errmin
;  fakei =
;  errtable.zmag_auto-errtable.imagxpsf+errtable.zmagxpsf-24.86663+25.67849
  fakei = x-c-24.86663+25.67849
  ierr = p[0]+p[1]*(fakei-20)+p[2]*(fakei-20)^2 > errmin
  w=where(fakei le 20)
  if w[0] ne -1 then ierr[w] = errmin
  ierr *= sqrt(obj->sxpar('EXPTIME', band='z')/obj->sxpar('EXPTIME', band='i'))
  out = sqrt(zerr^2+ierr^2)
  return, out
end
