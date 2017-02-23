Function find_rs, obj, _extra=extra
  cmd = cmdclass(obj, _extra=extra)
  wfit = where(cmd.rad and ~cmd.star and cmd.morph and cmd.rsregion)
  errs = geterrs(obj, _extra=extra)
  return, fit_rs(cmd.m[wfit], cmd.c[wfit], errs[wfit])
end
