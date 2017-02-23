Function bc03sspspec, $
   age, $
   metal, $
   imf=imf, $
   res=res, $
   pickles=pickles, $
   full=full

;  on_error, 2
  common jem$_bc03sspspec, _bc03sss

  if n_elements(imf) eq 0 then imf='chab'
  if n_elements(res) eq 0 then res='hr'

  ;; Pickles requested
  if keyword_set(pickles) then begin
     if metal ne 0.02 then message, 'Pickles lib only for Z=Z_solar'
     if n_elements(_bc03sss) eq 0 $
        || ~member(strupcase('pickles_'+imf),tag_names(_bc03sss)) then begin
        ;; need to load pickles data
        bc03dir=getenv('JEM_REF')+'SED/BC03/'
        if (file_info(bc03dir+'/sav/pickles_'+imf+'.sav')).exists then begin
           restore, bc03dir+'/sav/pickles_'+imf+'.sav'
        endif else begin
           file = bc03dir+'bc2003_hr_m62_'+imf+'_ssp_Pickles_Stelib.fits.gz'
           a=mrdfits(file,1,/silent)
           save, a, filename=bc03dir+'/sav/pickles_'+imf+'.sav'
        endelse
        if n_elements(_bc03sss) eq 0 then $
           _bc03sss = create_struct( 'pickles_'+imf, ptr_new(a) ) $
        else $
           _bc03sss = create_struct( _bc03sss, 'pickles_'+imf, ptr_new(a) )
     endif
     case imf of
        'chab' : ptr = _bc03sss.pickles_chab
        'salp' : ptr = _bc03sss.pickles_salp
     endcase
     tabinv, (*ptr).age, age, ageindex
     flux = interpolate( (*ptr).flux, ageindex )
     out = {wave:(*ptr).wave, flux:flux}
     if keyword_set(full) then begin
        tags = tag_names(*ptr)
        for i=0, n_elements(tags)-1 do begin
           if member(tags[i], ['WAVE','FLUX']) then continue
           out = create_struct(out, tags[i], $
                               interpolate( (*ptr).(i), ageindex ))
        endfor
     endif
  endif else begin

     ;; non-Pickles requested
     if n_elements(_bc03sss) eq 0 $
        || ~member(strupcase(imf+'_'+res),tag_names(_bc03sss)) then begin
        bc03dir=getenv('JEM_REF')+'SED/BC03/'
        if (file_info(bc03dir+'sav/'+imf+'_'+res+'.sav')).exists then begin
           restore, bc03dir+'sav/'+imf+'_'+res+'.sav'
        endif else begin
           files = file_search(bc03dir+'*_'+res+'*_'+imf+'*ssp.fits.gz', $
                               count=nfiles)
           a=mrdfits(files[0],1,/silent)
           ;; handle flux first as a special case and create structure aa
           nwave=n_elements(a.wave)
           nage=n_elements(a.age)
           aa=create_struct('AGE', a.age, $
                            'WAVE', a.wave, $
                            'FLUX', rebin(a.flux, nwave, nage, nfiles))
           tags = tag_names(a)
           for j=0, n_elements(tags)-1 do begin
              counter, j+1, n_elements(tags), 'Creating tag '
              if member(tags[j], ['AGE','WAVE','FLUX']) then continue
              aa=create_struct(aa,tags[j],rebin(a.(j), nage, nfiles))
           endfor
           for i=1, nfiles-1 do begin
              counter, i+1, nfiles
              b=mrdfits(files[i],1,/silent)
              tags = tag_names(a)
              for j=0, n_elements(tags)-1 do begin
                 if tags[j] eq 'FLUX' then begin
                    aa.flux[*,*,i] = b.flux
                    continue
                 endif
                 if member(tags[j], ['AGE','WAVE']) then continue
                 aa.(j)[*,i] = b.(j)
              endfor
           endfor
           save, aa, filename=bc03dir+'/sav/'+imf+'_'+res+'.sav'
        endelse
        if n_elements(_bc03sss) eq 0 then $
           _bc03sss = create_struct( imf+'_'+res, ptr_new(aa) ) $
        else $
           _bc03sss = create_struct( _bc03sss, imf+'_'+res, ptr_new(aa) )
     endif
     case imf of
        'chab' : begin
           case res of
              'hr' : ptr = _bc03sss.chab_hr
              'lr' : ptr = _bc03sss.chab_lr
           endcase
        end
        'salp' : begin
           case res of
              'hr' : ptr = _bc03sss.salp_hr
              'lr' : ptr = _bc03sss.salp_lr
           endcase
        end
     endcase
     metals = [0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05]
     tabinv, (*ptr).age, age, ageindex
     tabinv, metals, metal, metalindex
     flux = interpolate((*ptr).flux, $
                        ageindex, metalindex, $
                        /grid)
     out = {wave:(*ptr).wave, flux:flux}
     if keyword_set(full) then begin
        tags = tag_names(*ptr)
        for i=0, n_elements(tags)-1 do begin
           if member(tags[i], ['AGE','WAVE','FLUX']) then continue
           out = create_struct(out, tags[i], $
                               interpolate( (*ptr).(i), $
                                            ageindex, metalindex, $
                                            /grid ))
        endfor
     endif
  endelse
  return, out
end

Function bc03sspages, $
   imf=imf, $
   res=res, $
   pickles=pickles

   common jem$_bc03sspspec

  if n_elements(imf) eq 0 then imf='chab'
  if n_elements(res) eq 0 then res='hr'
  pickles = keyword_set(pickles)

  junk = bc03sspspec( 1.e9, 0.02, imf=imf, res=res, pickles=pickles )
  if keyword_set(pickles) then begin
     case imf of
        'chab' : ptr = _bc03sss.pickles_chab
        'salp' : ptr = _bc03sss.pickles_salp
     endcase
  endif else begin
     case imf of
        'chab' : begin
           case res of
              'hr' : ptr = _bc03sss.chab_hr
              'lr' : ptr = _bc03sss.chab_lr
           endcase
        end
        'salp' : begin
           case res of
              'hr' : ptr = _bc03sss.salp_hr
              'lr' : ptr = _bc03sss.salp_lr
           endcase
        end
     endcase
  endelse
  ages = (*ptr).age
  return, ages
end

Function sspages, $
   lib=lib, $
   _extra=extra

  if n_elements(lib) eq 0 then lib='BC03'

  case lib of
     'BC03' : ages = bc03sspages( _extra=extra )
  endcase

  return, ages
end

;; Load an SSP library from BC03 for now.  Will add Pegase, CB07,
;; Maraston in the future.  Possibly DustEM too?
Function sspspec, $
   age, $
   metal, $
   lib=lib, $
   _extra=extra

  if n_elements(lib) eq 0 then lib='BC03'

  case lib of
     'BC03' : spec = bc03sspspec( age, metal, _extra=extra )
  endcase

  return, spec
end
