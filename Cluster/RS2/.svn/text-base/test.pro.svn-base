Pro test
  print, 'Loading RDCS1252-2927 color-magnitude information'
  restore, '/home/scpdata02/cluster3/Y/Y.sav'
  resolve_obj, obj
  clustercmd = cmd2( obj, /xpsf )
  selectgals, clustercmd
  rMpc = 0.65
  rarcmin = rMpc / (lumdist( 1.237, /silent )/(1.+1.237)^2) * 180./!pi * 60.
  clusterarea = rarcmin^2*!pi
  mstar = 22.7
  
  print, 'Loading GOODS color-magnitude information'
  restore, '/home/scpdata03/goods/cmd/cmds.sav'
  bkggals = *(*cmds[0]).gals
  for i=1, 29 do begin
     cmd1 = *cmds[i]
     bkggals = [bkggals, *cmd1.gals]
  endfor
  bkgcmd = {clusterid:'G', clustername:'GOODS', zcluster:0.d, $
            fit:ptr_new(), gals:ptr_new(bkggals)}
  selectgals, bkgcmd
  bkgarea = 1.4^2*!pi*30.
  
  CMRparams = {intercept23:0.9, slope:-0.03, sigma:0.03}
  bkgparams = {d:2.0, e:0., f:0.}
  N = 30.

  print, rslike( clustercmd, bkgcmd, $
                 clusterarea, bkgarea, $
                 mstar+[-2.2,0.8], [0.6,1.2], $
                 CMRparams, bkgparams, N )
  
end
