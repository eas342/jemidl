Function TabularFactory::Init
  return, 1
end

Function TabularFactory::CreateFromFile, $
   filename
  filetype = self->inferFileType( filename )
  case filetype of 
     'fits': data = self->TabularFactory::CreateFromFitsFile( filename )
     'ascii': data = self->TabularFactory::CreateFromAsciiFile( filename )
  endcase
  return, data
end

Function TabularFactory::InferFileType, $
   filename
  if strpos( filename, 'fits' ) ne -1 then filetype = 'fits' $
  else if strpos( filename, 'dat' ) ne -1 then filetype = 'ascii' $
  else if strpos( filename, 'txt' ) ne -1 then filetype = 'ascii' $
  else filetype = ''
  return, filetype
end

Function TabularFactory::CreateFromAsciiFile, $
   filename
  ;;assume that the first line that doesn't start with '#' will set the format...
  a='#'
  openr, lun, filename, /get_lun
  while a[0] eq '#' do begin
     readf, lun, a
  endwhile
  aa = strsplit( a, /extract )
  ncolumns = n_elements( aa )
  
  case ncolumns of 
     1: begin
        close, lun
        free_lun, lun
        readcol, filename, format='D', column1
        data = column1
     end
     2: begin
        close, lun
        free_lun, lun
        readcol, filename, format='D,D', column1, column2
        data = [[column1],[column2]]
     end
     3: begin
        close, lun
        free_lun, lun
        readcol, filename, format='D,D,D', column1, column2, column3
        data = [[column1],[column2],[column3]]
     end
     else: begin
        data = double( aa )
        while not eof( lun ) do begin
           readf, lun, a
           if a[0] eq '#' then continue
           data = [[data], [double( strsplit( a, /extract ) )]]
        endwhile
        data = transpose(data)
        close, lun
        free_lun, lun
     end
  endcase
  return, transpose(data)
end

Function TabularFactory::CreateFromFitsFile, $
   filename, $
   columnwave=columnwave, $
   column1=column1, $
   column2=column2

  if n_elements( columnwave ) eq 0 then columnwave = ['WAVELENGTH', 'WAVE']
  if n_elements( column2 ) eq 0 then column2 = ['STATERROR', 'IVAR']
  if n_elements( column1 ) eq 0 then column1 = ['FLUX', 'THROUGHPUT']

  h=headfits( filename )
  if total( strpos( h, 'CRVAL1' ) ne -1 ) eq 0 then begin ;;no CRVAL1 in header
     a=mrdfits( filename, 1, /silent )
     tags = tag_names( a )
     for icw=0, n_elements( columnwave ) - 1 do begin
        w = (where( columnwave[icw] eq tags ))[0]
        if w ne -1 then begin
           wave = a.(w)
           break
        endif
     endfor
     if icw eq n_elements( columnwave ) then columnwave = '' $
     else columnwave = columnwave[icw]
     for ic1=0, n_elements( column1 ) - 1 do begin
        w = (where( column1[ic1] eq tags ))[0]
        if w ne -1 then begin
           c1data = a.(w)
           break
        endif
     endfor
     if ic1 eq n_elements( column1 ) then column1 = '' $
     else column1 = column1[ic1]
     for ic2=0, n_elements( column2 ) - 1 do begin
        w = (where( column2[ic2] eq tags ))[0]
        if w ne -1 then begin
           c2data = a.(w)
           break
        endif
     endfor
     if ic2 eq n_elements( column2 ) then column2 = '' $
     else column2 = column2[ic2]
  endif else begin
     ;;CRVAL1 data is available...
     naxis1 = sxpar( h, 'NAXIS1' )
     cdelt1 = sxpar( h, 'CDELT1' )
     cd1_1  = sxpar( h, 'CD1_1' )
     if cdelt1 eq 0. then cdelt1 = cd1_1
     crval1 = sxpar( h, 'CRVAL1' )
     wave = crval1 + dindgen( naxis1 )*cd1_1
     c1data = mrdfits( filename, 0, /silent )
     c2data = mrdfits( filename, 1, /silent )
  endelse
  data = [[wave],[c1data]]
  if n_elements( c1data ) eq n_elements( c2data ) then $
     data = [data, [c2data]]
  return, transpose(data)
end

Pro TabularFactory__Define
  struct = { TabularFactory, $
           a:0 }
end
