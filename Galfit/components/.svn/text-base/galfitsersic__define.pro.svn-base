Function GalfitSersic::Init, $
   _extra=extra
  
  self->SetProperty, fitxpos=1, fitypos=1, fitmag=1, fitre=1, fitn=1, fitboa=1, fitpa=1
  self->SetProperty, _extra=extra
  return, 1
end

Pro GalfitSersic::WriteToControlFile, lun, objnum
  printf, lun, string( format='(%"# Object number: %4i")', objnum )
  printf, lun, ' 0) sersic                       #  object type'
  printf, lun, string( format='(%" 1) %10.4f%10.4f%3i%3i   #  position x, y")', $
                       self.xpos, self.ypos, self.fitxpos, self.fitypos )
  printf, lun, string( format='(%" 3)%9.4f%3i                  #  integrated magnitude")', $
                       self.mag, self.fitmag )
  printf, lun, string( format='(%" 4)%9.4f%3i                  #  re (half-light radius)   [pix]")', $
                       self.re, self.fitre )
  printf, lun, string( format='(%" 5)%9.4f%3i                  #  Sersic index n (de Vaucouleurs n=4")', $
                       self.n, self.fitn )
  printf, lun, string( format='(%" 9)%9.4f%3i                  #  axis ratio (b/a)")', $
                       self.boa, self.fitboa )
  printf, lun, string( format='(%"10)%11.4f%3i                #  position angle (PA) [deg: Up=0, Left=90]")', $
                       self.pa, self.fitpa )
  printf, lun, string( format='(%" Z) %3i                          #  output option (0 = resid., 1 = Dont subtract)")', $
                       self.out )
  printf, lun
;  printf, lun, '# Object number: '+string(objnum, format='(I4)')
;  printf, lun, ' 1) '+string(self.xpos, format='(F10.4)')+string(self.ypos, format='(F10.4)')+string(self.fitxpos, format='(I3)')+string(self.fitypos, format='(I3)')+'   #  position x, y'
;  printf, lun, ' 3) '+string(self.mag, format='(F8.4)')+string(self.fitmag, format='(I3)')+'                  #  integrated magnitude'
;  printf, lun, ' 4) '+string(self.re,
;  format='(F8.4)')+string(self.fitre, format='(I3)')+'
;  #  re (half-light radius)   [pix]'
;  printf, lun, ' 5) '+string(self.n,
;  format='(F8.4)')+string(self.fitn, format='(I3)')+'
;  #  Sersic index n (de Vaucouleurs n=4'
;  printf, lun, ' 9) '+string(self.boa,
;  format='(F8.4)')+string(self.fitboa, format='(I3)')+'
;  #  axis ratio (b/a)'
;  printf, lun, '10) '+string(self.pa,
;  format='(F10.4)')+string(self.fitpa, format='(I3)')+'
;  #  position angle (PA) [deg: Up=0, Left=90]'
;  printf, lun, ' Z) '+string(self.out, format='(I3)')+"                          #  output option (0 = resid., 1 = Don't subtract)"
end

Function GalfitSersic::Results
  struct = { xpos: self.xpos, $
             ypos: self.ypos, $
             mag: self.mag, $
             re: self.re, $
             n: self.n, $
             boa: self.boa, $
             pa: self.pa, $
             errxpos: self.errxpos, $
             errypos: self.errypos, $
             errmag: self.errmag, $
             errre: self.errre, $
             errn: self.errn, $
             errboa: self.errboa, $
             errpa: self.errpa }
  return, struct
end

Pro GalfitSersic__Define, struct
  struct = { GalfitSersic, $
             Inherits GalfitComponent, $
;            Inherits JEMobject, $ ;; inherits via GalfitComponent
             xpos:0.d, $
             ypos:0.d, $
             mag:1.d, $
             re:1.d, $
             n:1.d, $
             boa:1.d, $
             pa:0.d, $
             fitxpos:1, $
             fitypos:1, $
             fitmag:1, $
             fitre:1, $
             fitn:1, $
             fitboa:1, $
             fitpa:1, $
             out:0, $
             errxpos:0.d, $
             errypos:0.d, $
             errmag:0.d, $
             errre:0.d, $
             errn:0.d, $
             errboa:0.d, $
             errpa:0.d }
end
