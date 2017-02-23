Pro testrslike
  restore, '/home/scpdata02/cluster3/R/R_2.sav'
  resolve_obj, obj
  s=obj->summary()
  clustercmd = cmd( obj, /xpsf, /nofit )
  clustergals = *clustercmd.gals

  restore, '/home/scpdata03/goods/GOODS-N-N01/N-N01.sav'
  obj->UpdateSpecCat
  bkgcmd = cmd( obj, /xpsf, /nofit )
  bkggals = *bkgcmd.gals

  clusterdata = {mag:0., color:0., msigma:0., csigma:0.}
  clusterdata = replicate(clusterdata, n_elements(clustergals))
  clusterdata.mag = clustergals.zmag
  clusterdata.color = clustergals.iz
  clusterdata.msigma = 0.
  clusterdata.csigma = clustergals.iz_err
  wcluster = where(clustergals.rad and ~clustergals.star)
  clusterdata = clusterdata[wcluster]

  bkgdata = {mag:0., color:0., msigma:0., csigma:0.}
  bkgdata = replicate(bkgdata, n_elements(bkggals))
  bkgdata.mag = bkggals.zmag
  bkgdata.color = bkggals.iz
  bkgdata.msigma = 0.
  bkgdata.csigma = bkggals.iz_err
  wbkg = where(bkggals.rad and ~bkggals.star)
  bkgdata = bkgdata[wbkg]

  clusterarea = 2.56
  bkgarea = 1.69
  magrange = [19, 24]
  colorrange = [0.7,1.1]
  CMRparams = { sigma:0.05, slope:-0.05, intercept23:1.0 }
  bkgparams = { d:0.1, e:0., f:0., g:0. }
  Schechterparams = { phi_star:1., m_star:22., alpha:-1. }
  iskysigma = mean(s.isky_sigma, /nan)
  zskysigma = mean(s.zsky_sigma, /nan)
  print, rslike( clusterdata, bkgdata, $
                 clusterarea, bkgarea, $
                 magrange, colorrange, $
                 iskysigma, zskysigma, $
                 CMRparams, bkgparams, Schechterparams )
end
