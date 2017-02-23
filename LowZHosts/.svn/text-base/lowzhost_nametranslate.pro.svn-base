Function LowZHost_nametranslate, name
  case name of
     'Ellip2' : return, 'SNF20080614-010'
     'Ellip3' : return, 'SNF20060521-008'
     'Ellip4' : return, 'SNF20070417-002'
     'Ellip5' : return, 'SNF20070727-016'
     'Ellip6' : return, 'SNF20080731-000'
     'PTF10dnp' : return, 'PTF10fps'
     else:
  endcase
  if strmid(name, 0, 2) eq '20' then return, 'SN'+name
  if strmid(name, 0, 2) eq 'sn' then return, 'SN'+strmid(name, 2)
  return, name
end
