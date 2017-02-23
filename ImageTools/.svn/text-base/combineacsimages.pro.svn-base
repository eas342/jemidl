;+
; NAME:
;	COMBINEACSIMAGES
;
; PURPOSE:
;       This procedure takes a list of single exposure ACS images and
;       attempts to coadd them.
; CALLING SEQUENCE:
;	COMBINEACSIMAGES, image_list, outfilename
; INPUTS:
;       image_list: A string vector of single exposure ACS filenames
;       outfilename: Where to place the output.
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       The file [outfilename] will be created.
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       JEM, Jan, 2009. Written
;-

Pro CombineACSimages, image_list, filename
  nimages=n_elements(image_list)
  ;; need exposure times and image size
  exptime=fltarr(nimages)
  for i=0, nimages-1 do begin
     h=headfits(image_list[i])
     if i eq 0 then begin
        nx=sxpar(h,'NAXIS1')
        ny=sxpar(h,'NAXIS2')
     endif
     exptime[i]=sxpar(h,'EXPTIME')
  endfor
  ;; create output vars
  outputimage=fltarr(nx,ny)
  ;; for now use 4x4 grid processing
  npartitions=1
  xpartitions=intarr(npartitions,2)
  ypartitions=intarr(npartitions,2)
  ;; calculate partition begin and end pixels
  for i=0,npartitions-2 do begin
     x=floor(nx*((i+1.)/npartitions))
     y=floor(ny*((i+1.)/npartitions))
     xpartitions[i,1]   = x
     xpartitions[i+1,0] = x+1
     ypartitions[i,1]   = y
     ypartitions[i+1,0] = y+1
  endfor
  xpartitions[npartitions-1,1]=nx-1
  ypartitions[npartitions-1,1]=ny-1

; LOOP OVER PARTITIONS
  for xpart=0, npartitions-1 do begin
     for ypart=0, npartitions-1 do begin
        ;; recall partition begin/end
        xsize=xpartitions[xpart,1]-xpartitions[xpart,0]+1
        ysize=ypartitions[ypart,1]-ypartitions[ypart,0]+1
        ;; temp storage for subimages
        images=fltarr(xsize,ysize,nimages)
        time=fltarr(xsize,ysize,nimages)
        errimage=fltarr(xsize,ysize,nimages)
        ;; loop over images
        for iimage=0, nimages-1 do begin
        ;;for iimage=0, 9 do begin
           print, xpart, ypart, iimage
           a=mrdfits(image_list[iimage],/silent)
           h=headfits(image_list[iimage])
           gain=mean( double( [sxpar( h, 'ATODGNA' ), $
                               sxpar( h, 'ATODGNB' ), $
                               sxpar( h, 'ATODGNC' ), $
                               sxpar( h, 'ATODGND' )] ) )
           rdnoise = mean( double( [sxpar( h, 'READNSEA' ), $
                                    sxpar( h, 'READNSEB' ), $
                                    sxpar( h, 'READNSEC' ), $
                                    sxpar( h, 'READNSED' )] ) )
           sky = double( sxpar( h, 'SKY_MEAN' ) )  ;; counts per 500s(?)
           sky /= 500                              ;; seconds per exposure
           
           dark=0.00318 ;; counts per second
           
           ;; only keep relevant subimages
           images[*,*,iimage]=a[xpartitions[xpart,0]:xpartitions[xpart,1], $
                                ypartitions[ypart,0]:ypartitions[ypart,1]]
           ;; and exptime subimage
           a=mrdfits(strmid(image_list[iimage],0,strlen(image_list[iimage])-8)+'wht.fits',/silent)
           time[*,*,iimage]=a[xpartitions[xpart,0]:xpartitions[xpart,1], $
                                ypartitions[ypart,0]:ypartitions[ypart,1]]
           tmp = time[*,*,iimage] lt 0.05*max(time[*,*,iimage])
           for ishift=-1, 1 do begin
              for jshift=-1, 1 do begin
                 time[*,*,iimage] *= shift(1-tmp, ishift, jshift)
              endfor
           endfor
           errimage[*,*,iimage] = ErrIm(images[*,*,iimage],time[*,*,iimage], gain, rdnoise, sky, dark, 1)
        endfor
        ;; now do the averaging
        print, 'averaging!'
        w=where(~finite(errimage))
        if w[0] ne -1 then errimage[w]=0.
        a=jem_avsigclip(images, 1./errimage^2, 2., 2., 10)
        outputimage[xpartitions[xpart,0]:xpartitions[xpart,1], $
                    ypartitions[ypart,0]:ypartitions[ypart,1]] $
           = a
     endfor
  endfor
;; output results
  mwrfits, outputimage, filename, /create
end
