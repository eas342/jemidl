Function helio_corr, hdr
  telescope = strcompress(strmid(sxpar(hdr[*, 0], 'TELESCOP'), 0, 3) $
                          , /rem)
  if telescope eq '' $
     and strcompress(sxpar(hdr[*,0], 'INSTRUME'),/rem) eq 'OSIRIS' then telescope='GTC'

  case telescope OF
     'Kec': begin  ;; Keck/LRIS
        mjd = double(sxpar(hdr, 'MJD-OBS'))  + 2400000.5D
        equinox = double(sxpar(hdr, 'EQUINOX'))
        ra = sxpar(hdr, 'RA')
        dec = sxpar(hdr, 'DEC')
        x_radec, ra, dec, radeg, decdeg
     end
     'GTC': begin ;; GTC
        mjd = sxpar(hdr, 'MJD-OBS')+2400000.5D
        equinox = 2000.0
        obs = 'lapalma'
        ra = sxpar(hdr, 'RA')
        dec = sxpar(hdr, 'DEC')
        x_radec, ra, dec, radeg, decdeg
     end
  endcase
  helio = -1.0*x_keckhelio(radeg, decdeg, equinox, jd = mjd, OBS = OBS)
  return, helio
end
