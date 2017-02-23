Pro testclean2
  restore, '/home/jmeyers314/olivetos/jemidl/CLEAN/Apsf.sav'
  a=mrdfits('/home/jmeyers314/scp2/clusters2/R/galfit/R1154z.fits',1)
  b=clean(a.image,zpsf,1.5)
  atv, [a.image,b]
end
