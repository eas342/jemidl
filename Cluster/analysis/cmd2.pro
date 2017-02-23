Function cmd2, $
   obj, $
   _extra=extra

  resolve_obj, obj
  s=obj->summary()
  iz=obj->color( _extra=extra, /silent )
  z850=s.zmag_auto
  id=s[0].clusterid
  clustername=obj->extract('clustername')
  zcluster=obj->extract('zcluster')

  ;; MW dust correction
  ebv = obj->extract('ebv')
  A_i = 1.973*ebv
  A_z = 1.472*ebv
  z850 -= A_z
  iz -= (A_i - A_z)

  ;; SExtractor MAG_AUTO bias correction
  z850 += z850*(-0.0422781)+0.635273

  izerr=geterrs7( obj, ierr=ierr, zerr=zerr )

  ;; Find clustercentric radius
  ra_center=obj->extract('ra')
  dec_center=obj->extract('dec')
  ra=s.alphawin_j2000
  dec=s.deltawin_j2000
  radius2 = (ra_center-ra)^2*(cos(dec_center*!dpi/180.d))^2 $
            + (dec_center-dec)^2 ;; clustercentric distance^2 in degrees^2
  rad = sqrt(radius2)*!dpi/180.d ;; in radians
  rad *= lumdist( zcluster, /silent )/(1.+zcluster)^2 ;; in Mpc

  ;; Put results in structure
  @define_structs2
  gals = replicate( CMDgal0, n_elements(s) )
  gals.galid = s.galid
  gals.starcheck = s.zflux_radius lt 2.2
  gals.no2 = s.o2_ew le -10.
  gals.ra = s.alphawin_j2000
  gals.dec = s.deltawin_j2000
  gals.aim = s.a_image
  gals.bim = s.b_image
  gals.kron = s.kron_radius
  gals.theta = s.theta_image
  gals.rad = rad
  gals.z850 = z850
  gals.z850err = zerr
  gals.i775 = z850 + iz 
  gals.i775err = ierr
  gals.iz = iz
  gals.izerr = izerr
  gals.gini = s.gini
  gals.asym = s.asym
  gals.re = s.re
  gals.n = s.n
  gals.zspec = s.z
  gals.o2ew = s.o2_ew
  gals.o2ewerr = s.o2_ewerr


  ;; Identify SNe hosts in structure
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

  ;; Create CMD structure
  cmd = cmdcluster0
  cmd.clusterid = id
  cmd.clustername = clustername
  cmd.zcluster = zcluster
  cmd.gals = ptr_new(gals)

  return, cmd
end
