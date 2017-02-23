Pro HTMLheader, lun, id
  printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN"'
  printf, lun, '    "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">'
  printf, lun, '<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" >'
  printf, lun, '<head>'
  printf, lun, '<meta http-equiv="Content-Type" content="application/xhtml+xml; charset=UTF-8" />'
  printf, lun, '<title> Josh Meyers | Cluster-'+id+' CMD</title>'
  printf, lun, '<meta name="Author" content="Josh Meyers" />'
  printf, lun, '<style type="text/css">'
  printf, lun, 'body {font-family: tahoma; arial, sans-serif;}'
  printf, lun, '#imap {display:block; width:1000px; height:500px; background:url(CMD.jpg) no-repeat; position:relative;}'
  printf, lun, '#imap a#title {display:block; width:1000px; height:0; padding-top:500px; overflow:hidden; position:absolute; left:0; top:0; background:transparent url(CMD.jpg) no-repeat 1000px 1000px; cursor:default;}'
  printf, lun, '#imap a#title:hover {background-position: 0 0; z-index:10;}'
end

Pro MakeWebCMD, obj, dir, _extra=extra
  dir = directoryify(dir)
  file_mkdir, dir
  resolve_obj, obj
  window, 0, xsize=1000, ysize=500, /pixmap
;  plot1cmd1, obj, /specmem, /specnon, /sne, /clean
  plot_cmd, obj, /specmem, /specnon, /sne, _extra=extra
  image=tvread()
  write_jpeg, dir+'CMD.jpg', image, true=1
  
  s=obj->summary()
  mags = s.zmag_auto
  colors = obj->color( _extra=extra )
  ebv = obj->extract('ebv')
  A_i = 1.973*ebv
  A_z = 1.472*ebv
  mags -= A_z
  colors -= (A_i - A_z)
  mags += mags*(-0.0422781)+0.635273

  openw, lun, dir+'index.html', /get_lun
  id = obj->extract('clusterid')
  HTMLheader, lun, id

  
  ra=obj->extract('ra')
  dec=obj->extract('dec')
  z=obj->extract('zcluster')
  radius = sqrt(((ra-s.alphawin_j2000)^2*(cos(dec*!dpi/180.d))^2 $
                 + (dec-s.deltawin_j2000)^2))
  radius *= (lumdist( z, /silent )/(1.+z)^2)*!dpi/180.d
  cmd = cmdclass(obj, _extra=extra)
  wgals = where( mags gt 19 and mags lt 25 and colors gt -0.5 and colors lt 1.5 $
                 and ~cmd.star and cmd.rad)
  printf, lun, '#imap dd {position:absolute; padding:0; margin:0;}'
  for igal=0, n_elements(wgals)-1 do begin
     pix = convert_coord(mags[wgals[igal]], colors[wgals[igal]], /data, /to_device)
     left = string(pix[0]-4, format='(I4)')
     top = string(500-(pix[1]+5), format='(I4)')
     name = id+string(wgals[igal],format='(I05)')
     printf, lun, '#imap #'+name+' {left:'+left+'px; top:'+top+'px; z-index:20;}'
  endfor

  for igal=0, n_elements(wgals)-1 do begin
     name = id+string(wgals[igal],format='(I05)')
     if igal eq 0 then str = '#imap a#'+name+'id' $
     else str = strjoin([str, ', #imap a#'+name+'id'])
  endfor
  printf, lun, strjoin([str, ' {display:block; width:10px; height:10px; text-decoration:none; z-index:20;}'])

  printf, lun, '#imap a span, #imap a:visited span {display:none;}'

  for igal=0, n_elements(wgals)-1 do begin
     name = id+string(wgals[igal],format='(I05)')
     if igal eq 0 then str = '#imap a#'+name+'id:hover' $
     else str = strjoin([str, ', #imap a#'+name+'id:hover'])
  endfor
  printf, lun, strjoin([str, ' {background-position:0 0}'])

  printf, lun, '#imap a:hover span {position:absolute;  width:300px; display:block; font-family:arial; font-size:12px; background:#fff; color:#000; border:1px solid #000; padding:5px;}'
  printf, lun, '#imap a span:first-line { font-weight:bold; font-style:italic; }'

  
  for igal=0, n_elements(wgals)-1 do begin
     pix = convert_coord(mags[wgals[igal]], colors[wgals[igal]], /data, /to_device)
     name = id+string(wgals[igal],format='(I05)')
     left = string(-(pix[0]-4), format='(I4)')
     top = string(pix[1]+8, format='(I4)')
     printf, lun, '#imap a#'+name+'id:hover span {left:'+left+'px; top:'+top+'px;}'
  endfor

  printf, lun, '</style>'
  printf, lun, '</head>'

  printf, lun, '<body>'
  printf, lun, '<dl id="imap">'
  printf, lun, ' <dt><a id="title" href="#nogo" title="CMD">CMD</a></dt>'

  obj->LoadImage
  for igal=0, n_elements(wgals)-1 do begin
     counter, igal+1, n_elements(wgals)
     s1 = s[wgals[igal]]
     name = id+string(wgals[igal],format='(I05)')
     printf, lun, ' <dd id="'+name+'"><a id="'+name+'id" title="'+name+'" href="#nogo"><span>'
     printf, lun, '  '+name
     printf, lun, string(format='(%"  <br />z<sub>850</sub>: %05.2f       i<sub>775</sub> - z<sub>850</sub>: %+05.2f")', mags[wgals[igal]], colors[wgals[igal]])
     printf, lun, string(format='(%"  <br />RA: %10.6f    DEC: %+09.5f")', s1.alphawin_j2000, s1.deltawin_j2000)
     printf, lun, string(format='(%"  <br />asym: %6.3f  gini: %6.3f")', s1.asym2, s1.gini2)
;     printf, lun, string(format='(%"  <br />asym: %6.3f  conc: %6.3f  gini: %6.3f")', s1.asym, s1.conc, s1.gini)
     printf, lun, string(format='(%"  <br />R: %6.3fMpc")', radius[wgals[igal]])
;     printf, lun, string(format='(%"  <br />fluxrad: %6.3fpx")', s1.zflux_radius)
     if s1.z ge 0. then $
        printf, lun, string(format='(%"  <br />z: %6.3f")', s1.z)
     fileprefix = string(format='(%"%s%05i")', id, wgals[igal])

     igalimage = obj->galimage(wgals[igal], band='i', xsize=50, ysize=50)
     write_png, dir+fileprefix+'i.png', bytscl(igalimage, min=-0.01, max=0.1)
     printf, lun, '  <br /><img src="'+fileprefix+'i.png'+'">'

     zgalimage = obj->galimage(wgals[igal], band='z', xsize=50, ysize=50)
     write_png, dir+fileprefix+'z.png', bytscl(zgalimage, min=-0.01, max=0.1)
     printf, lun, '  <img src="'+fileprefix+'z.png'+'">'

     image = uindgen(4,51,51)
     image[0,*,*] = 256*bytscl(zgalimage, min=-0.01, max=0.09)
     image[1,*,*] = 256*bytscl(0.5*(zgalimage+igalimage), min=-0.005, max=0.1)
     image[2,*,*] = 256*bytscl(igalimage, min=-0.005, max=0.1)
     image[3,*,*] = 65535L
     write_png, dir+fileprefix+'c.png', image
     printf, lun, '  <img src="'+fileprefix+'c.png'+'">'
     printf, lun, ' <br /></span></a></dd>'
  endfor
  obj->FreeImages
  printf, lun, '</dl>'
  printf, lun, '</body>'
  printf, lun, '</html>'
  close, lun
  free_lun, lun
  wdelete, 0
end
