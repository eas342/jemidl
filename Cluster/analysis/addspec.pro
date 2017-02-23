Pro addspec, obj
  resolve_routine, 'spectrum__define'
  resolve_routine, 'dopplershift__define'
  resolve_obj, obj
  id = obj->extract('clusterid')
  case id of 
     'A' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Hilton2009/', 'Aspec.cat'
        ;;obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-A/'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/VLT/CL-A/', 'Aspec.cat'
     end
     'B' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-B/', 'Bspec.cat'
     end
     'C' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/VLT/CL-C/', 'Cspec.cat'
     end
     'D' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-O/', 'Ospec.cat'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/VLT/CL-O/', 'Ospec.cat'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Andreon2008/', 'Dspec.cat'
     end
     'E' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-E/', 'Espec.cat'
     end
     'F' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-F/', 'Fspec.cat'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-F/', 'Fspec.cat'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Eisenhardt2008/', 'Fspec.cat'
     end
     'G' : begin
        ;;obj->addspectroscopy,'/home/jmeyers314/scp1/CL-spec/Subaru/CL-G'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Eisenhardt2008/', 'Gspec.cat'
     end
     'H' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-H/HA/', 'HAspec.cat'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-H/HB/', 'HBspec.cat'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-H/', 'Hspec.cat'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Brodwin2006/', 'Hspec.cat'
     end
     'I' : begin
        ;;obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-I/'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Eisenhardt2008/', 'Ispec.cat'
     end
     'J' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Eisenhardt2008/', 'Jspec.cat'
     end
     'K' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-K/', 'Kspec.cat'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-K/', 'Kspec.cat'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Stanford2005/', 'Kspec.cat'
     end
     'L' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-L/', 'Lspec.cat'
     end
     'M' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Postman1998/', 'Mspec.cat'
     end
     'N' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-N/', 'Nspec.cat'
     end
     'P' : begin
     end
     'R' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-R/RA/', 'RAspec.cat'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-R/RB/', 'RBspec.cat'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-R/RC/', 'RCspec.cat'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/VLT/CL-R/', 'Rspec.cat'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Bremer2006/', 'Rspec.cat'
     end
     'S' : begin
     end
     'T' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-T/', 'Tspec.cat'
     end
     'U' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/VLT/CL-U/', 'Uspec.cat'
     end
     'V' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Keck/CL-V/', 'Vspec.cat'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-V/', 'Vspec.cat'
     end
     'W' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Rosati1999/', 'Wspec.cat'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Private/StanfordW/', 'Wspec.cat'
     end
     'X' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Subaru/CL-X/', 'Xspec.cat'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Stanford2002/', 'Xspec.cat'
     end
     'Y' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/DeMarco2007/', 'Yspec.cat'
     end
     'Z' : begin
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/VLT/CL-Z/', 'Zspec.cat'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Rosati_b/', 'Zspec1.cat'
        obj->addspectroscopy, '/home/jmeyers314/scp1/CL-spec/Literature/Rosati_b/', 'Zspec2.cat'
     end
     else :
  endcase
  obj->UpdateSpecCat
end
