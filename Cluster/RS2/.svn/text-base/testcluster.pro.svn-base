pro testcluster, id
  readcol, '/home/scpdata02/cluster3/cmd2/fit_rs7/CMD.txt', $
           ids, mstar, int23, slope, ngal, $
           format='A,F,F,X,X,F,X,X,X,X,F'
  w=where(id eq ids)
  slope=slope[w[0]]
  mstar=mstar[w[0]]
  int23=int23[w[0]]
  ngal=ngal[w[0]]
  restore, '/home/scpdata02/cluster3/'+id+'/'+id+'.sav'
  resolve_obj, obj
  start = [int23, 0.07, ngal, 1.5, 0.]
  result = rsmetropoliswrapper( obj, slope, magrange=mstar+[-2.2,0.8], $
                                rsresidrange=int23+[-0.3,0.3], /plothist, $
                                start=start )
  
end
