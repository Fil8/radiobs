#!/usr/bin/env python

import sys, os, string
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table

import fluxInt

fint=fluxInt.flInt()

# Print help message and exit
if 'help' in sys.argv or '-h' in sys.argv:
    print 'Run as follows'
    print 'fluxint.py -f <filename.fits> OR -t <filename.tbl> AND -m <maskname.fits> OR -r <region.fits>'
    print ' OPTIONS (-c <cutoff> -fr <observed_frequency in MHz> -flerrr <systematic error on flux (%)>)'
    sys.exit()
else: arg=sys.argv

if '-f' in arg:
    fileName=arg[arg.index('-f')+1]
elif '-t' in arg: 
    tableName=arg[arg.index('-t')+1]
    tableFileNames= ascii.read(tableName,format='csv')

if '-m' in arg:
    maskName=arg[arg.index('-m')+1]
    maskData,maskHead = fint.openFile(maskName)
    maskData = np.squeeze(maskData)
    maskHead = fint.cleanHead(maskHead)
    fits.writeto(maskName,maskData,maskHead,overwrite=True)
elif '-r' in arg:
    region =arg[arg.index('-r')+1]

if '-c' in arg:
    cutoff =float(arg[arg.index('-c')+1])
else:
    print 'WARNING: using default cutoff 1e-3*units of image'
    cutoff = 1e-3

if '-flerr' in arg:
    errFlux = float(arg[arg.index('-flerr')+1])
else:
    print 'WARNING: using default error on Flux: 0.05*integrated flux'
    errFlux = 5


if ('-f' in arg and '-r' in arg):
    print '\t Executing File+Region Combo' 
    datas,heads = fint.openFile(fileName)
    datas=np.squeeze(datas)
    heads=fint.cleanHead(heads)
    maskedData, background, noise, pixels=fint.maskDatReg(datas,heads,region,cutoff)

    freq=fint.readFreq(heads,arg)
    fluxint, number_pix_beam=fint.measFlux(maskedData,heads,noise,errFlux)
    t,flErr = fint.writeTable(heads,fluxint,noise,number_pix_beam,freq,errFlux)

    print t

if ('-f' in arg and '-m' in arg):
    print '\t Executing File+Mask Combo' 
    maskedData, background, noise, pixels, newHeads = fint.maskData(fileName,maskName)
    heads=newHeads

    freq=fint.readFreq(heads,arg)
    
    fluxint, numPixBeam=fint.measFlux(maskedData,heads,noise,errFlux)
    
    t,flErr = fint.writeTable(heads,fluxint,noise,numPixBeam,freq,errFlux)

    print t

if ('-t' in arg and '-m' in arg):
    print '\t Executing FileList+Mask Combo'
    for i in xrange(0,len(tableFileNames.columns[0])):

        fileName = tableFileNames.columns[0][i]
        maskedData, background, noise, pixels, newHeads = fint.maskData(fileName,maskName)
        heads=newHeads
        
        if 'Frequency' in tableFileNames.dtype.names:
            freq = tableFileNames.columns[1][i]*1e6

        fluxint, numPixBeam =fint.measFlux(maskedData,heads,noise,errFlux)

        t,flErr = fint.writeTable(heads,fluxint,noise,numPixBeam,freq,errFlux)
        print t

if ('-t' in arg and '-r' in arg):
    print '\t Executing FileList+Region Combo' 
    for i in xrange(0,len(tableFileNames.columns[0])):

        fileName = tableFileNames.columns[0][i]
        datas,heads = fint.openFile(fileName)
        datas=np.squeeze(datas)
        heads=fint.cleanHead(heads)
        maskedData, background, noise, pixels=fint.maskDatReg(datas,heads,region,cutoff)
        
        if 'Frequency' in tableFileNames.dtype.names:
            freq = tableFileNames.columns[1][i]*1e6

        fluxint, numPixBeam =fint.measFlux(maskedData,heads,noise,errFlux)

        t,flErr = fint.writeTable(heads,fluxint,noise,numPixBeam,freq,errFlux)

        print t