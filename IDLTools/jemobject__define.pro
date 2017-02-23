Pro JEMobject::SetProperty, _extra=extra
  if n_elements(extra) eq 0 then return
  properties=tag_names(extra)

  call_procedure, obj_class(self)+'__define', struct
  for i=0L, n_tags(extra)-1 do begin
     theProperty = properties[i]
     index = where(strpos(tag_names(struct), theproperty) eq 0, count)
     index = index[0]
     if count gt 1 then message, "Ambiguous keyword: " $
                                 +theProperty $
                                 +". Use more characters in it's specification."
     if count eq 0 then message, "Keyword not found."
     self.(index) = extra.(i)
  endfor
end

Function JEMobject::Extract, field
  s=size(field)
  if s[s[0]+1] ne 7 then message, 'Field variable must be a string'
  thisClass = obj_class(self)
  ok = execute('thisStruct = {'+thisClass+'}')
  structFields = Tag_Names(thisStruct)
  index = where(structFields eq StrUpCase(field), count)
  if count eq 1 then begin
     if (execute('retVal = self.'+structFields[index[0]])) then begin
        if n_elements(retVal) gt 1 then return, retVal
        if (ptr_valid(retVal)) then begin
           return, *retVal
        endif else begin
           return, retVal
        endelse
     endif
     return, -1
  endif else begin
     message, 'Can not find field "'+field+ $
              '" in object '+obj_class(self), /informational
     return, -1
  endelse
end

Pro JEMobject__Define, struct
  struct = { JEMobject, $
             junk:0L }
end
