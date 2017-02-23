Function GalfitSky::Init, _extra=extra
  self->SetProperty, fitsky=0, fitdskydx=0, fitdskydy=0
  self->SetProperty, _extra=extra
  return, 1
end

Pro GalfitSky::WriteToControlFile, lun, objnum
  printf, lun, '# Object number: '+string(objnum, format='(I4)')
  printf, lun, ' 0) sky                          # object type'
  printf, lun, ' 1) '+string(self.sky, format='(E17.4E2)')+string(self.fitsky, format='(I3)')+'                  # sky background at center of fitting region [ADUs]'
  printf, lun, ' 2) '+string(self.dskydx, format='(F9.5)')+string(self.fitdskydx, format='(I3)')+'                  # dsky/dx (sky gradient in x)'
  printf, lun, ' 3) '+string(self.dskydy, format='(F9.5)')+string(self.fitdskydy, format='(I3)')+'                  # dsky/dy (sky gradient in y)'
  printf, lun, ' Z) '+string(self.out, format='(I3)')+"                          #  output option (0 = resid., 1 = Don't subtract)"

end

Function GalfitSky::Results
  struct = { sky:       self.sky, $
             dskydx:    self.dskydx, $
             dskydy:    self.dskydy, $
             errsky:    self.errsky, $
             errdskydx: self.errdskydx, $
             errdskydy: self.errdskydy }
  return, struct
end

Pro GalfitSky__Define, struct
  struct = { GalfitSky, $
             Inherits GalfitComponent, $
;             Inherits JEMobject, $  ;; inherits via GalfitComponent
             sky:       0.d, $
             dskydx:    0.d, $
             dskydy:    0.d, $
             fitsky:    0, $
             fitdskydx: 0, $
             fitdskydy: 0, $
             out:       0, $
             errsky:    0.d, $
             errdskydx: 0.d, $
             errdskydy: 0.d }
end
