Function GalfitAbsConstraint::Init, $
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

Pro GalfitAbsConstraint::WriteToConstraintFile, lun
  printf, lun, format='(%"%i  %s  %f to %f")', self.objnum, self.param, self.min, self.max
end

Pro GalfitAbsConstraint__define
  struct = { GalfitAbsConstraint, $
             Inherits GalfitConstraint, $
             objnum:0L, $
             param:'', $
             min:0.d, $
             max:0.d }
end
