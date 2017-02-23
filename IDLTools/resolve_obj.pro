;+
; NAME:
;	RESOLVE_OBJ
;
; PURPOSE:
;       This routine will restore class methods for objects
;       restored with save.
; CALLING SEQUENCE:
;	RESOLVE_OBJ, object
; INPUTS:
;       object: An object restored from an IDL save file for which
;       you wish to restore class methods.
; KEYWORD PARAMETERS:
;       None.
; OUTPUTS:
;       None.
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       Taken from http://www.dfanning.com/tips/saved_objects.html
;-

Pro Resolve_Obj, obj, class=class, routine_info=ri
  if n_params() ne 0 then begin
     if not obj_valid(obj) then begin
        message, 'Object is not valid.'
     endif
     class=obj_class(obj)
  endif

  if n_elements(ri) eq 0 then ri=routine_info()

  for i=0, n_elements(class)-1 do begin
     defpro=class[i]+'__DEFINE'
     if (where(ri eq defpro))[0] eq -1 then begin
        call_procedure, defpro
     endif
     supers=obj_class(class[i],/superclass, count=cnt)
     if cnt gt 0 then resolve_obj, class=supers, routine_info=ri
  endfor
end
