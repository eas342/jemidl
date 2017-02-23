Pro testy
  restore, '/home/scpdata02/cluster3/Y/Y_2.sav'
  resolve_obj, obj
  clustercmd = cmd( obj, /xpsf, /nofit )
  clustergals = *clustercmd.gals
  s=obj->summary()
  
  restore, '/home/scpdata03/goods/N-N01.sav'
  obj->UpdateSpecCat
  bkgcmd = cmd( obj, /xpsf, /nofit )
  bkggals = *bkgcmd.gals

  colorrange = [0.6, 1.2]
  magrange = [20, 24.]

  clusterarea = 2.56
  bkgarea = 1.69
  CMRparams = { sigma:0.05, slope:-0.05, intercept23:1.0 }
  bkgparams = { d:0.1, e:0., f:0., g:0. }
  Schechterparams = { phi_star:1., m_star:22., alpha:-1. }
  iskysigma = mean(s.isky_sigma, /nan)
  zskysigma = mean(s.zsky_sigma, /nan)


end
