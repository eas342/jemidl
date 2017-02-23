Function imagesection, image, p, background=background
  type=size( image, /tname )
  case type of
     'BYTE'     : backdef = 0
     'INT'      : backdef = 0
     'UINT'     : backdef = 0
     'LONG'     : backdef = 0
     'ULONG'    : backdef = 0
     'LONG64'   : backdef = 0
     'ULONG64'  : backdef = 0
     'COMPLEX'  : backdef = complex( !values.f_nan, !values.f_nan )
     'DCOMPLEX' : backdef = dcomplex( !values.d_nan, !values.d_nan )
     'FLOAT'    : backdef = !values.f_nan
     'DOUBLE'   : backdef = !values.d_nan
  endcase
  if n_elements(background) eq 0 then background=backdef
  sz=size( image, /dim )
  out = fltarr(p[2]-p[0]+1, p[3]-p[1]+1)+background
  if p[0] ge sz[0] or p[1] ge sz[1] or p[2] lt 0 or p[3] lt 0 then return, out
  p_out = [0, 0, p[2]-p[0], p[3]-p[1]]
  p_in = p
  if p[0] lt 0 then begin
     p_in[0] = 0
     p_out[0] = -p[0]
  endif
  if p[1] lt 0 then begin
     p_in[1] = 0
     p_out[1] = -p[1]
  endif
  if p[2] ge sz[0] then begin
     p_in[2] = sz[0]-1
     p_out[2] = sz[0]-p[0]-1
  endif
  if p[3] ge sz[1] then begin
     p_in[3] = sz[1]-1
     p_out[3] = sz[1]-p[1]-1
  endif
  out[p_out[0]:p_out[2], p_out[1]:p_out[3]] = image[p_in[0]:p_in[2], p_in[1]:p_in[3]]
  return, out
end
