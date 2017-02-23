Function struct_renametags, struct, newtags
  struct1 = struct[0]
  ntags = n_tags(struct)
  oldtags = tag_names(struct)
  newstruct = create_struct(newtags[ntags-1], struct1.(ntags-1))
  for itag=1, ntags-1 do begin
     newstruct = create_struct(newtags[ntags-itag-1], struct1.(ntags-itag-1), newstruct)
  endfor
  dims = size(struct,/dim)
  ndim=n_elements(dims)
  case ndim of
     1: newstruct = replicate(newstruct, dims[0])
     2: newstruct = replicate(newstruct, dims[0], dims[1])
     3: newstruct = replicate(newstruct, dims[0], dims[1], dims[2])
  endcase
  for itag=0, ntags-1 do begin
     newstruct.(itag) = struct.(itag)
  endfor
  return, newstruct
end
