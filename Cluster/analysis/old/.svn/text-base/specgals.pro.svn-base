Function specgals, objs
  resolve_obj, objs[0]
  for iobj=0, n_elements(objs)-1 do begin
     speccat=objs[iobj]->Extract('speccat')
     if size(speccat,/tname) ne 'STRUCT' then continue
     sexcat=objs[iobj]->Extract('sexcat')
     galapagoscat=objs[iobj]->Extract('galapagoscat')
     galfitcat=objs[iobj]->Extract('galfitcat')
     morphcat=objs[iobj]->Extract('morphcat')
     
     z = objs[iobj]->extract('zcluster')
     veldisp = objs[iobj]->extract('veldisp')
     if z eq 0. then continue
     ckms=299792.d
     struct = { zcluster:z, $
                veldisp:veldisp, $
                velocity:0.d }
     struct = replicate(struct, n_elements(sexcat))
     struct.velocity=(speccat.z-z)*ckms/(1.+z)
     outcat1 = struct_addtags(sexcat, struct)
     outcat1 = struct_addtags(outcat1, struct_trimtags(galapagoscat, except=['CLUSTERID','GALID']))
     outcat1 = struct_addtags(outcat1, struct_trimtags(galfitcat, except=['CLUSTERID','GALID']))
     outcat1 = struct_addtags(outcat1, struct_trimtags(speccat, except=['CLUSTERID','GALID']))
     outcat1 = struct_addtags(outcat1, struct_trimtags(morphcat, except=['CLUSTERID','GALID']))

     w = where( speccat.zqual ne 'N' $
                and speccat.zqual ne 'F' $
                and speccat.zqual ne 'C' $
                and sexcat.zclass_star le 0.7 $
                and sexcat.ellipticity le 0.5 $
                and sexcat.zfwhm_image ge 3.5 )
     if n_elements(outcat) eq 0 then outcat = outcat1[w] $
     else outcat=[outcat, outcat1[w]]
  endfor
  return, outcat
end
