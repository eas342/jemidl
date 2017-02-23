Function FilterFactory::Init
  if not self->TabularFactory::Init() then return, 0
  return, 1
end

Function FilterFactory::CreateFromFile, $
   filename, $
   johnson=johnson

  data = self->TabularFactory::CreateFromFile( filename )
  if keyword_set(johnson) then begin
     data[1,*] /= data[0,*]
     data[1,*] /= max(data[1,*])
  endif
  filter = obj_new('filter', data[0,*], data[1,*] )
  return, filter
end

Function FilterFactory::CreateFromFitsFile, $
   filename, $
   johnson=johnson

  data = self->TabularFactory::CreateFromFitsFile( filename )
  if keyword_set(johnson) then begin
     data[1,*] /= data[0,*]
     data[1,*] /= max(data[1,*])
  endif
  filter = obj_new('filter', data[0,*], data[1,*] )
  return, filter
end

Function FilterFactory::CreateFromAsciiFile, $
   filename, $
   johnson=johnson

  data = self->TabularFactory::CreateFromAsciiFile( filename )
  if keyword_set(johnson) then begin
     data[1,*] /= data[0,*]
     data[1,*] /= max(data[1,*])
  endif
  filter = obj_new('filter', data[0,*], data[1,*] )
  return, filter
end

Pro FilterFactory__Define
  struct = { FilterFactory, $
             Inherits TabularFactory, $
             c:0 }
end
