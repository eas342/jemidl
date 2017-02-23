;; execute this in various Final/ and FinalCR/ subdirectories
pro doindiv_index_check
  cd, current=pwd
  case 1 of
     (strpos(pwd, 'blue600_4000_d460') ne -1) : begin
        setup='B460'
     end
     (strpos(pwd, 'blue600_4000_d560') ne -1) : begin
        setup='B560'
     end
     (strpos(pwd, 'red600_5000') ne -1) : begin
        setup='R5000'
     end
     (strpos(pwd, 'red600_7500') ne -1) : begin
        setup='R7500'
     end
     (strpos(pwd, 'R1000B') ne -1) : begin
        setup='R1000B'
     end
     (strpos(pwd, 'R2000B') ne -1) : begin
        setup='R2000B'
     end
     (strpos(pwd, 'R2500R') ne -1) : begin
        setup='R2500R'
     end
  endcase
  f=file_search('*check.fits')
  for i=0, n_elements(f)-1 do begin
     indiv_index_check, f[i], setup
  endfor
end
