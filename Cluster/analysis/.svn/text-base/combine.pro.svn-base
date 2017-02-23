Pro combine, rootdir, objs=objs
  
  readcol, 'clusters.txt', ID, name, name2, ra, dec, nz, zcluster, vel_disp, format='A,A,A,A,A,I,D,D'
  dirs=file_search(rootdir+'/?',/mark_dir)
  junk=obj_new('cluster')

  for idir=0, n_elements(dirs)-1 do begin
;  for idir=0, 0 do begin
     id = strmid(dirs[idir],strlen(dirs[idir])-2,1)
     filename = dirs[idir]+id+'.sav'
     if ~(file_info(filename)).exists then continue
     print, 'restoring '+filename
     restore, filename, /relax
     resolve_obj, obj

     case id of 
        'A' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Hilton/'
           ;;obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-A/'
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/VLT/CL-A/'
        end
        'B' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-B/'
        end
        'C' : begin
           ;;obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-C/'
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/VLT/CL-C/'
        end
        'D' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-O/'
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/VLT/CL-O/'
        end
        'E' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-E/'
        end
        'F' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-F/'
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-F/'
        end
        'G' : begin
           ;;obj->addspectroscopy,'/home/jmeyers314/scp1/CL-spec/Subaru/CL-G'
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Eisenhardt_G'
        end
        'H' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-H/HA/'
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-H/HB/'
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-H/'
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Brodwin/'
        end
        'I' : begin
           ;;obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-I/'
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Eisenhardt_I'
        end
        'J' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Eisenhardt_J'
        end
        'K' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-K/'
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-K/'
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Stanford_b/'
        end
        'L' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-L/'
        end
        'M' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Postman/'
        end
        'N' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-N/'
        end
        'P' : begin
           ;;obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-P/'
        end
        'R' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-R/RA/'
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-R/RB/'
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-R/RC/'
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/VLT/CL-R/'           
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Bremer/'
        end
        'S' : begin
           ;;obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-S/'
        end
        'T' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-T/'
        end
        'U' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/VLT/CL-U/'
        end
        'V' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-V/'
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-V/'
        end
        'W' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Rosati/'
        end
        'X' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-X/'
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Stanford/'
        end
        'Y' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/DeMarco/'
        end
        'Z' : begin
           obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/VLT/CL-Z/'
        end
        else :
     endcase
     if n_elements(objs) eq 0 then objs=[obj] $
     else objs=[objs,obj]
     obj->SpecRegionCheck, 'reg/'+id+'spec.reg'
;     obj->AddMorphology, dirs[idir]+'galfit/'
  endfor
end
