Pro testclean
  restore, 'Apsf.sav'
  a=mrdfits('/home/jmeyers314/scp2/clustergalaxies/R/galfit/R1154z.fits',1)
  clean=obj_new('CLEAN',a.image, zpsf, 1.86, loopgain=0.5)
  thresh=6.
  val = 0.003
  for i=0, 3000 do begin
     void=clean->Clean(1000, absval=val)
     clean->print
     a=clean->stats()
     print, [a,a[0]/a[1]]
     if a[0] lt val then break
  endfor
  print, i*1000
  obj_destroy, clean
end
