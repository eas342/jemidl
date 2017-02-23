Function GalfitRelConstraint::Init, $
   objnum, $
   param, $
   min, $
   max
  self.objnum = objnum
  self.param = param
  self.min = min
  self.max = max
  return, 1
end

Pro GalfitRelConstraint::WriteToConstraintFile, lun
  printf, lun, format='(%"%i  %s  %f %f")', self.objnum, self.param, self.min, self.max
end

Pro GalfitRelConstraint__define
  struct = { GalfitRelConstraint, $
             Inherits GalfitConstraint, $
             objnum:0L, $
             param:'', $
             min:0.d, $
             max:0.d }
end
