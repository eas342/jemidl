Pro RSreg
  rootdir='/home/scpdata03/cluster5/'
  readcol, rootdir+'catalogs/clusters.txt', ids, name1, name2, ra, dec, z, $
           format='A,A,A,A,A,F', /silent
  for i=0, n_elements(ids)-1 do begin
     id=ids[i]
     catdir = rootdir+'catalogs/CL-'+id+'/'
     s = mrdfits( catdir+'CL-'+id+'_summary.fits', 1, /silent )
     wRS = where( s.gini gt 0.4 and s.asym lt 0.1 and s.ccradmpc lt 0.65 and abs(s.clrsresidfixm[2]/s.clrsresidfixmerr[2]) lt 2 )

     openw, lun, catdir+'CL-'+id+'.RS.reg', /get_lun
     printf, lun, 'global color=green font="helvetica 10 normal" '+ $
             'select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source'
     printf, lun, 'fk5'
     printf, lun
     printf, lun, format='(%"circle( %s, %s, %s )")', ra[i], dec[i], mpc2arcmin( 0.65, z[i] )/60.
     if wRS[0] ne -1 then begin
        for j=0, n_elements(wRS)-1 do begin
           if s[wRS[j]].clrsresidfixm[2] gt 0 then color='red' else color='magenta'
           printf, lun, format='(%"circle( %s, %s, %s ) # color=%s text={%7.3f}")', $
                   s[wRS[j]].alpha_j2000, s[wRS[j]].delta_j2000, s[wRS[j]].g2re/3600./20., color, s[wRS[j]].iz
        endfor
     endif
     close, lun
     free_lun, lun
  endfor
end
