Function like_cluster, data, CMRparams
  sigpart = CMRparams.sigma^2 + data.izerr^2
  constpart = 1./sqrt(2.*!pi*sigpart)
  exppart = -(data.rsresid - CMRparams.intercept23)^2
  exppart /= 2.*sigpart
  return, constpart * exp(exppart)
end

Function like_field, data, bkgparams
  return, bkgparams.d $
          + bkgparams.e*data.rsresid $
          + bkgparams.f*data.rsresid^2
end

Function bkgintegrand, rsresid
  common JEM$_rslike4, bkgparams, CMRparams, p, q
  return, like_field({rsresid:rsresid}, bkgparams)
end

Function clusterintegrand, rsresid
  common JEM$_rslike4
  return, like_cluster({rsresid:rsresid, izerr:0.}, $
                       CMRparams )
end

Function rslike, clustercmd, bkgcmd, $
                 clusterarea, bkgarea, $
                 magrange, rsresidrange, $
                 CMRparams1, bkgparams1, N
  common JEM$_rslike4
  CMRparams=CMRparams1
  bkgparams=bkgparams1

  clustergals = *clustercmd.gals
  clusterrsresid = clustergals.iz - CMRparams.slope*(clustergals.z850-23.)
  wcluster =  where( clustergals.select $
                     and clustergals.z850 ge magrange[0] $
                     and clustergals.z850 le magrange[1] $
                     and clusterrsresid ge rsresidrange[0] $
                     and clusterrsresid le rsresidrange[1] )
  clusterdata = { rsresid:clusterrsresid[wcluster], $
                  izerr:clustergals[wcluster].izerr }
  
  bkggals = *bkgcmd.gals
  bkgrsresid = bkggals.iz - CMRparams.slope*(bkggals.z850-23.)
  wbkg = where( bkggals.select $
                and bkggals.z850 ge magrange[0] $
                and bkggals.z850 le magrange[1] $
                and bkgrsresid ge rsresidrange[0] $
                and bkgrsresid le rsresidrange[1] )
  bkgdata = { rsresid:bkgrsresid[wbkg], $
              izerr:bkggals[wbkg].izerr }
  
  clusterpart = like_cluster( clusterdata, CMRparams )
  fieldpart = like_field( clusterdata, bkgparams )
  bkgpart = like_field( bkgdata, bkgparams )
  
  bkgintegral = qromb( 'bkgintegrand', $
                       rsresidrange[0], rsresidrange[1] )
  clusterintegral = N*qromb( 'clusterintegrand', $
                           rsresidrange[0], rsresidrange[1] )

  print, clusterintegral
  return, total(alog(clusterpart*N+fieldpart)) + alog(clusterarea) $
          +total(alog(bkgpart))+alog(bkgarea) $
          - bkgintegral*(clusterarea+bkgarea) $
          - clusterintegral*clusterarea
end
