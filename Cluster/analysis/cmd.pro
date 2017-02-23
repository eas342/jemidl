Function cmd, obj, colorrange=colorrange, magrange=magrange, $
              colorfitrange=colorfitrange, bluelimit=bluelimit, $
              zrange=zrange, slope=slope, morph=morph, $
              clusterrad=clusterrad, plotrange=plotrange, $
              nofit=nofit, err=err, nsigma=nsigma, $
              nfitsigma=nfitsigma, intercept=intercept, $
              iscatter=iscatter, fixslope=fixslope, $
              pslope=pslope, pwid=pwid, fit5=fit5, $
              _extra=extra
  
  if n_elements(clusterrad) eq 0 then clusterrad = 0.8 ;;Mpc
  if n_elements(morph) eq 0 then morph=8 ;;morphology criteria #8
  if n_elements(err) eq 0 then err=7
  if n_elements(slope) eq 0 then slope=0.
  if n_elements(bluelimit) eq 0 then bluelimit=0.
  if n_elements(colorrange) eq 0 then colorrange=[0.8, 1.1]
  if n_elements(magrange) eq 0 then magrange=[20.,23.5]
  if n_elements(colorfitrange) eq 0 then colorfitrange=colorrange
  if n_elements(iscatter) eq 0 then iscatter=0.03
  if ~finite(iscatter) then iscatter = 0.03

  resolve_obj, obj
  s=obj->summary()
  iz=obj->color( _extra=extra, /silent )
  zmag_auto=s.zmag_auto
  id=obj->extract('clusterid')

  ;; MW dust correction
  ebv = obj->extract('ebv')
  A_i = 1.973*ebv
  A_z = 1.472*ebv
  zmag_auto -= A_z
  iz -= (A_i - A_z)
  
  ;; SExtractor MAG_AUTO bias correction
  ;m -= (m-21)*0.1/3.+0.2
  zmag_auto += zmag_auto*(-0.0422781)+0.635273

  ;; clustercentric radius cut
  zcluster=obj->extract('zcluster')
  ra_center=obj->extract('ra')
  dec_center=obj->extract('dec')
  ra=s.alphawin_j2000
  dec=s.deltawin_j2000
  radius2 = (ra_center-ra)^2*(cos(dec_center*!dpi/180.d))^2 $
            + (dec_center-dec)^2 ;; clustercentric distance^2 in degrees^2
  rthresh = clusterrad/(lumdist( zcluster, /silent )/(1.+zcluster)^2) ;; in rad
  rthresh *= 180.d/!dpi ;; rad->degrees
  radcheck = radius2 le rthresh^2

  ;; star cut
  starcheck = s.zflux_radius lt 2.2
  
  ;; morphcut
  case morph of
     1: morphcheck = s.conc ge 0.35+4.65*s.asym
     2: morphcheck = s.asym le 0.1 and s.gini ge 0.44
     3: morphcheck = s.asym le 0.08 and s.gini ge 0.44
     4: morphcheck = s.asym le 0.1 and s.gini ge 0.43
     5: morphcheck = s.asym le 0.1 and s.gini ge 0.4
     6: morphcheck = s.asym le 0.1 and s.gini ge 0.45
     7: morphcheck = s.asym le 0.13 and s.gini ge 0.4
     8: morphcheck = s.asym le 0.1 and s.gini ge 0.42
  endcase

  if n_elements(zrange) eq 0 then begin
     if strmid(obj->extract('clustername'), 0, 5) eq 'GOODS' then begin
        zrange = [0., 1.]
     endif else begin
        readcol, '/home/scpdata02/clusters/zlimits.txt', ids, zmin, zmax, format='A,F,F', /silent
        iid = where(ids eq id)
        zrange = [zmin[iid], zmax[iid]]
     endelse
  endif
  specmemcheck = s.z ge zrange[0] and s.z le zrange[1]
  specnoncheck = s.z ge 0. and ~specmemcheck
  no2check = strpos( s.comment, 'oII' ) eq -1 and (specmemcheck or specnoncheck) $
             and s.comment ne '?' and s.comment ne '??'
  o2check = strpos( s.comment, 'oII' ) ne -1
  nocluecheck = s.comment eq '?' or s.comment eq '??'
  case err of 
     3: iz_err=geterrs3(obj)
     4: iz_err=geterrs4(obj)
     5: iz_err=geterrs5(obj)
     6: iz_err=geterrs6(obj)
     7: iz_err=geterrs7(obj, ierr=ierr, zerr=zerr)
  endcase
  zmag_auto_err = s.zmagerr_auto

  if n_elements(nfitsigma) eq 0 then begin
     rsregioncheck = rsregioncheck( zmag_auto, iz, magrange, colorrange, bluelimit, slope=slope )
     rsfitcheck = rsregioncheck( zmag_auto, iz, magrange, colorfitrange, bluelimit, slope=slope )
  endif else begin
     rsregioncheck = rsnsigmacheck( zmag_auto, iz, iz_err, nsigma, magrange, bluelimit, $
                                    intercept=intercept, slope=slope, iscatter=iscatter )
     rsfitcheck = rsnsigmacheck( zmag_auto, iz, iz_err, nfitsigma, magrange, bluelimit, $
                                 intercept=intercept, slope=slope, iscatter=iscatter )
  endelse

  @define_structs
  gals = replicate(CMDgal0, n_elements(s))
  gals.galid = s.galid
  gals.rad = radcheck
  gals.star = starcheck
  gals.morph = morphcheck
  gals.rsregion = rsregioncheck
  gals.specmem = specmemcheck
  gals.specnon = specnoncheck
  gals.no2 = no2check
  gals.o2 = o2check
  gals.noclue = nocluecheck
  gals.zmag_auto = zmag_auto
  gals.zmag_auto_err = zmag_auto_err
  gals.iz = iz
  gals.iz_err = iz_err
  gals.zmag = zmag_auto
  gals.zmagerr = zerr
  gals.imag = zmag_auto+iz
  gals.imagerr = ierr
  gals.gini = s.gini
  gals.asym = s.asym
  gals.re = s.re
  gals.n = s.n
  gals.zspec = s.z

  readcol, '/home/scpdata02/clusters/supernova.txt', $
           SNename, SNnickname, SNtype, SNhostra, SNhostdec, $
           format='A,A,X,A,X,X,A,A', /silent
  wSNe = where( id eq strmid( SNename, 5, 1 ) )
  if wSNe[0] ne -1 then begin
     for iw=0, n_elements(wSNe)-1 do begin
        if SNhostra[wSNe[iw]] eq '??' then continue ;;hostless SN
        get_coords, coords, instring=SNhostra[wSNe[iw]]+' '+SNhostdec[wSNe[iw]]
        nearby = obj->nearestobjs(coords[0]*15., coords[1])
        gals[nearby[0]].SNname = SNename[wSNe[iw]]   
        gals[nearby[0]].SNtype = SNtype[wSNe[iw]]   
     endfor
  endif

  out = cmdcluster0
  out.clusterid = s[0].clusterid
  out.clustername = obj->extract('clustername')
  out.zcluster = zcluster
  out.colorrange = colorrange
  out.magrange = magrange
  out.slope = slope
  out.gals = ptr_new(gals)
                
  if ~keyword_set( nofit ) then begin
     wfit = where( gals.rad and $
                   ~gals.star and $
                   gals.morph and $
                   gals.rsregion )
     fit = where( rsfitcheck[wfit] )
     
     if keyword_set(fit5) then begin
        out.fit = ptr_new(fit_rs5( gals[wfit].zmag_auto, gals[wfit].iz, gals[wfit].iz_err, $
                                   fixslope=fixslope, int23=intercept ))
     endif else if n_elements(fixslope) ne 0 then begin
        out.fit = ptr_new(fit_rs3( gals[wfit].zmag_auto, gals[wfit].iz, gals[wfit].iz_err, $
                                   fit=fit, slope=fixslope ))
     endif else if n_elements(pslope) ne 0 then begin
        out.fit = ptr_new(fit_rs4( gals[wfit].zmag_auto, gals[wfit].iz, gals[wfit].iz_err, $
                                   pslope=pslope, pwid=pwid ))
     endif else begin 
        out.fit = ptr_new(fit_rs2( gals[wfit].zmag_auto, gals[wfit].iz, gals[wfit].iz_err, $
                                   fit=fit, iscatter=iscatter ))
     endelse
     
     w=where( gals.rad and ~gals.star and gals.morph)
     c=gals[w].iz-(*out.fit).slope*(gals[w].zmag_auto-25)
     plothist, c, xhist, yhist, peak=1, /noplot, xrange=[-0.5,1.5], bin=0.05, /nan
     out.xhist = ptr_new(xhist)
     out.yhist = ptr_new(yhist)
  endif

  return, out
end
