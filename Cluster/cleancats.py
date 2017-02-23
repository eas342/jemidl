def line2string(line) :
    str="%10i"%int(line[0])         # NUMBER
    str += "%11.3f"%float(line[1])  #XWIN_IMAGE
    str += "%11.3f"%float(line[2])  #YWIN_IMAGE
    str += "%12.7f"%float(line[3])  #ALPHA
    str += "%12.7f"%float(line[4])  #DELTA
    str += "%11.3f"%float(line[5])  #FLUX_RADIUS
    str += "%6.2f"%float(line[6])   #CLASS_STAR
    str += "%10.3f"%float(line[7])  #A_IMAGE
    str += "%10.3f"%float(line[8])  #B_IMAGE
    str += "%6.1f"%float(line[9])   #THETA_IMAGE
    str += "%9.3f"%float(line[10])  #ELLIPTICITY
    str += "%6.2f"%float(line[11])  #KRON_RADIUS
    str += "%6.2f"%float(line[12])  #PETRO_RADIUS
    str += "%9.4f"%float(line[13])  #MU_MAX
    str += "%9.4f"%float(line[14])  #MU_THRESHOLD
    str += "%14.10f"%float(line[15]) #BACKGROUND
    str += "%10i"%float(line[16])   #ISOAREA_IMAGE
    str += "%10.3f"%float(line[17]) #FWHM_IMAGE
    str += "%9.4f"%float(line[18])  #MAG_BEST
    str += "%9.4f"%float(line[19])  #MAGERR_BEST
    str += "%9.4f"%float(line[20])  #MAG_AUTO
    str += "%9.4f"%float(line[21])  #MAGERR_AUTO
    str += "%9.4f"%float(line[22])  #MAG_ISO
    str += "%9.4f"%float(line[23])  #MAGERR_ISO
    str += "%9.4f"%float(line[24])  #MAG_ISOCOR
    str += "%9.4f"%float(line[25])  #MAGERR_ISOCOR
    str += "%9.4f"%float(line[26])  #MAG_PETRO
    str += "%9.4f"%float(line[27])  #MAGERR_PETRO
    for l in line[28:] :            #MAG_APER + MAGERR_APER
        if float(l) > 99. : str += "%9.4f"%99.
        else : str += "%9.4f"%float(l)
    str += '\n'
    return str

import subprocess
import sys

if len(sys.argv) < 3 :
    print 'usage: '
    print 'python cleancats.py i.cat z.cat'
    sys.exit()

ifile = sys.argv[1]
zfile = sys.argv[2]

fzin = open(zfile)
fiin = open(ifile)

fzout = open(zfile+'.new','w')
fiout = open(ifile+'.new','w')

zlines = fzin.readlines()
ilines = fiin.readlines()

iobj = 0
print len(ilines)
for iline,zline in zip(*[ilines,zlines]) :
    if iline[0] == '#' : 
        fiout.write(iline)
        fzout.write(zline)
        continue
    istuff=iline.split()
    zstuff=zline.split()
    if istuff[19] == 'inf' : #19 = MAGERR_BEST
        continue
    if float(istuff[19]) > 3. :
        continue
    if float(zstuff[1]) < 0. or float(zstuff[2]) < 0. : #XWIN,YWIN.  This actually happened once in the Y cluster!
        continue
    if abs(float(istuff[1]) - float(zstuff[1])) > 8. : #cut out stuff that doesn't have a consistent center
        continue
    if abs(float(istuff[2]) - float(zstuff[2])) > 8. : #cut out stuff that doesn't have a consistent center
        continue
    if float(istuff[23]) > 40. : #23 = MAGERR_ISO
        istuff[23] = str(99.)
    if float(istuff[25]) > 40. : #25 = MAGERR_ISOCOR
        istuff[25] = str(99.)
    if float(istuff[27]) > 40. : #27 = MAGERR_PETRO
        istuff[27] = str(99.)
    istuff[0] = iobj
    zstuff[0] = iobj
    #print iline
    #print line2string(istuff)
    #print zline
    #print line2string(zstuff)
    fzout.write(line2string(zstuff))
    fiout.write(line2string(istuff))
    iobj += 1
    
fzout.close()
fzin.close()
fiout.close()
fiin.close()
cmd = "mv "+zfile+".new "+zfile
subprocess.call(cmd, shell=True)
cmd = "mv "+ifile+".new "+ifile
subprocess.call(cmd, shell=True)

