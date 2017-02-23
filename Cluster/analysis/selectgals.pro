Pro selectgals, cmd, $
                ginicut=ginicut, $
                asymcut=asymcut, $
                radcut=radcut, $
                zrange=zrange, $
                o2lim=o2lim
  if n_elements(ginicut) eq 0 then ginicut = 0.42
  if n_elements(asymcut) eq 0 then asymcut = 0.1
  if n_elements(radcut) eq 0 then radcut = 0.65 ;;Mpc
  if n_elements(zrange) eq 0 then begin
     if strmid(cmd.clustername, 0, 5) eq 'GOODS' then begin
        zrange = [0., 1.]
     endif else begin
        readcol, '/home/scpdata02/clusters/zlimits.txt', ids, zmin, zmax, format='A,F,F', /silent
        iid = where(ids eq cmd.clusterid)
        zrange = [zmin[iid], zmax[iid]]
     endelse
  endif
  if n_elements(o2lim) eq 0 then o2lim = -10.


  gals=*cmd.gals
  select = gals.gini ge ginicut $
           and gals.asym le asymcut $
           and gals.rad le radcut
  (*cmd.gals).select = select
  (*cmd.gals).morphcheck = gals.gini ge ginicut and gals.asym le asymcut
  (*cmd.gals).radcheck = gals.rad le radcut

  (*cmd.gals).specmem = gals.zspec ge zrange[0] and gals.zspec le zrange[1]
  (*cmd.gals).specnon = gals.zspec ge 0. and ~(*cmd.gals).specmem
  (*cmd.gals).o2 = gals.o2ew lt o2lim
  (*cmd.gals).no2 = gals.o2ew ge o2lim
  (*cmd.gals).noclue = gals.zspec ge 0. and ~finite(gals.o2ew)
  
end
