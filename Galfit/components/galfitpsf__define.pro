Function GalfitPSF::Init, $
   _extra=extra

  self->SetProperty, fitxpos=1, fitypos=1, fitmag=1
  self->SetProperty, _extra=extra
  return, 1
end

Pro GalfitPSF::WriteToControlFile, lun, objnum
  printf, lun, string( format='(%"# Object number: %4i")', objnum )
  printf, lun, ' 0) psf                          #  object type'
  printf, lun, string( format='(%" 1) %10.4f%10.4f%3i%3i   #  position x, y")', $
                       self.xpos, self.ypos, self.fitxpos, self.fitypos )
  printf, lun, string( format='(%" 3) %8.4f%3i                  #  integrated magnitude")', $
                       self.mag, self.fitmag )
  printf, lun, string( format='(%" Z) %3i                          #  output option (0 = resid., 1 = Dont subtract)")', self.out )
  printf, lun
end

Function GalfitPSF::Results
  struct = { xpos: self.xpos, $
             ypos: self.ypos, $
             mag:  self.mag, $
             errxpos: self.errxpos, $
             errypos: self.errypos, $
             errmag:  self.errmag }
  return, struct
end

Pro GalfitPSF__Define, struct
  struct = { GalfitPSF, $
             Inherits GalfitComponent, $
;            Inherits JEMobject, $ ;; inherits via GalfitComponent
             xpos:    0.d, $
             ypos:    0.d, $
             mag:     1.d, $
             fitxpos: 1, $
             fitypos: 1, $
             fitmag:  1, $
             out:     0, $
             errxpos: 0.d, $
             errypos: 0.d, $
             errmag:  0.d }
end
  
