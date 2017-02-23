Pro combine1, rootdir, objs=objs
  
  readcol, '/home/scpdata02/clusters/clusters.txt', ID, name, name2, ra, dec, nz, zcluster, vel_disp, format='A,A,A,A,A,I,F,F'
  dirs=file_search(rootdir+'/?',/mark_dir)
  junk=obj_new('cluster')

  for idir=0, n_elements(dirs)-1 do begin
     id = strmid(dirs[idir],strlen(dirs[idir])-2,1)
     filename = dirs[idir]+id+'.sav'
     if ~(file_info(filename)).exists then continue
     print, 'restoring '+filename
     restore, filename, /relax
     resolve_obj, obj
     addspec, obj
     obj->FreeImages
     obj->UpdateSpecCat
     if n_elements(objs) eq 0 then objs=[obj] $
     else objs=[objs,obj]
  endfor
end
