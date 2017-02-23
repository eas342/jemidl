Function SpectrumFactory::Init
  if not self->TabularFactory::Init() then return, 0
  return, 1
end

Function SpectrumFactory::CreateFromFile, $
   filename, $
   unit=unit, $
   orig_z=orig_z, $
   new_z=new_z, $
   _extra=extra


  if n_elements( unit ) eq 0 then unit = 'flambda'
  if n_elements( orig_z ) eq 0 then orig_z = obj_new( 'dopplershift', z=0.d )
  if n_elements( new_z ) eq 0 then new_z = obj_new( 'dopplershift', z=0.d )
  
  filetype = self->inferFileType( filename )
  case filetype of
     'fits': spectrum = self->CreateFromFitsFile( filename, unit=unit, $
                                                  orig_z=orig_z, new_z=new_z, $
                                                  _extra=extra )
     'ascii': spectrum = self->CreateFromAsciiFile( filename, unit=unit, $
                                                    orig_z=orig_z, new_z=new_z, $
                                                    _extra=extra )
  endcase
  return, spectrum
end

Function SpectrumFactory::CreateFromFitsFile, $
   filename, $
   unit=unit, $
   orig_z=orig_z, $
   new_z=new_z, $
   errfilename=errfilename, $
   varfilename=varfilename

  sdata = self->TabularFactory::CreateFromFitsFile( filename, column2=column2 )
  wave=sdata[0,*]
  flux=sdata[1,*]
  dims = size( sdata, /dimens )
  if dims[1] eq 3 then begin
     case column2 of 
        'STATERROR': ivar = 1.d/sdata[2,*]^2.d
        'IVAR': ivar = sdata[2,*]
     endcase
  endif
  
  if n_elements(errfilename) ne 0 then begin
     edata = self->TabularFactory::CreateFromFitsFile( errfilename, column2=['STATERROR'] )
     ivar = 1.d/edata[1,*]^2.d
  endif
  if n_elements(varfilename) ne 0 then begin
     vdata = self->TabularFactory::CreateFromFitsFile( varfilename, column2=['VAR','VARIANCE'] )
     ivar = 1.d/vdata[1,*]
  endif
  
  if n_elements(flux) ne n_elements(ivar) then $
     spectrum = obj_new('spectrum', wave, flux, unit=unit, orig_z=orig_z, new_z=new_z ) $
  else $
     spectrum = obj_new('spectrum', wave, flux, ivar=ivar, unit=unit, orig_z=orig_z, new_z=new_z )
  return, spectrum
end

Function SpectrumFactory::CreateFromAsciiFile, $
   filename, $
   errfilename=errfilename, $
   varfilename=varfilename, $
   unit=unit, $
   orig_z=orig_z, $
   new_z=new_z

  sdata = self->TabularFactory::CreateFromAsciiFile( filename )
  wave = sdata[0,*]
  flux = sdata[1,*]
  dims = size( sdata, /dimens )
  if dims[0] eq 3 then ivar = sdata[2,*]

  if n_elements(errfilename) ne 0 then begin
     edata = self->TabularFactory::CreateFromAsciiFile( errfilename )
     ivar = 1.d/edata[1,*]^2.d
  endif
  if n_elements(varfilename) ne 0 then begin
     vdata = self->TabularFactory::CreateFromAsciiFile( varfilename )
     ivar = 1.d/idata[1,*]
  endif

  if n_elements(flux) ne n_elements(ivar) then $
     spectrum = obj_new('spectrum', wave, flux, unit=unit, orig_z=orig_z, new_z=new_z ) $
  else $
     spectrum = obj_new('spectrum', wave, flux, ivar=ivar, unit=unit, orig_z=orig_z, new_z=new_z )
  return, spectrum
end

Pro SpectrumFactory__Define
  struct = { SpectrumFactory, $
             Inherits TabularFactory, $
             b:0 }
end
