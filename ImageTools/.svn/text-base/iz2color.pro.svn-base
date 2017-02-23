Pro iz2color, ifile, zfile, outfile
  blue=mrdfits(ifile)
  red=mrdfits(zfile, 0, hdr)
  green=0.5*(blue+red)
  sz=size(blue,/dim)
  out=fltarr( sz[0], sz[1], 3 )
  out[*,*,0]=red
  out[*,*,1]=green
  out[*,*,2]=blue
  mwrfits, out, outfile, hdr, /create
end
