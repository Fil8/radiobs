#!/usr/bin/env python

__author__ = "Filippo Maccagni"
__copyright__ = "Fil8"
__email__ = "filippo.maccagni@gmail.com"

import sys, os, string
import numpy as np
from astropy.io import fits
from astropy import units as u
from reproject import reproject_exact
from mpdaf.obj import Image
import fluxInt

fint = fluxInt.fluxint()

class cvRegr:
    

    def __init__(self):

        self.rootdir = os.getcwd()+'/'
        #self.filename = sys.argv[1]

    def convolve(self,basename,main_beam):
            
        outName = string.split(basename,'.fits')[0]
        outName = outName+'_cv.fits'            
        dat, baseheader = fint.openFile(basename)
        basefile = fits.open(basename)
        #baseheader = basefile[0].header
        basedata=Image(basename)
        
        if 'NAXIS3' in baseheader:
            del baseheader['NAXIS3']
        if 'NAXIS4' in baseheader:
            del baseheader['NAXIS4'] 
 
        aaa= basedata.data
        
        beam = np.array([float(baseheader['BMAJ']),float(baseheader['BMIN'])])
        print main_beam, beam
        if main_beam[0] > beam[0] and main_beam[1] > beam[1] :
            bx = np.sqrt(main_beam[0]*main_beam[0]-beam[0]*beam[0])
            by = np.sqrt(main_beam[1]*main_beam[1]-beam[1]*beam[1])
            #bx= main_beam[0]
            #by = main_beam[1]
            if 'CDELT1' in baseheader:
                pix_size = -baseheader['CDELT1']        
            elif 'CD1_1' in baseheader:       
                pix_size = -baseheader['CD1_1']
        
        
            beam_area = 2*np.pi*beam[0]/2.35482*beam[1]/2.35482
            main_beam_area = 2*np.pi*main_beam[0]/2.35482*main_beam[1]/2.35482

            number_pix_beam = beam_area/(pix_size*pix_size)

            #aaa = np.divide(aaa,number_pix_beam)
            #aaa = np.squeeze(aaa)  
            #basedata.data = np.squeeze(basedata.data)    
            #basedata.data=aaa
            #print basedata.shape

            newdata = Image.fftconvolve_gauss(basedata,center =None,flux=1, peak=True,factor=1,fwhm=(by,bx),
                    unit_center=u.degree, unit_fwhm=u.degree, inplace=False)    

            aaa = np.array(newdata.data)
        
            result = np.divide(aaa,np.power(main_beam_area/beam_area,2))


            #result = np.multiply(aaa,main_beam_area/(pix_size*pix_size))    
            
            baseheader['BMAJ'] = main_beam[0]
            baseheader['BMIN'] = main_beam[1]

            fits.writeto(outName,result,baseheader,overwrite=True)
        
        else:
            outName=basename

        return outName



    def regrid(self,basename,slavename):
        
        outName = string.split(slavename,'.fits')[0]
        outName = outName+'_rg.fits'            


        bdata, bheader = fint.openFile(basename)
        sdata, sheader = fint.openFile(slavename)
        slave = fits.open(slavename)[0]
        
        bheader['BMIN'] = sheader['BMIN']
        bheader['BMAJ'] = sheader['BMAJ']
        if 'FREQ' in slave.header:
            bheader['FREQ'] = sheader['FREQ']
        elif 'CRVAL3' in sheader:
            bheader['FREQ'] = sheader['CRVAL3']

        #print basename
        #for i in base.header.keys():
        #    print i,'\t',base.header[i]
        #print slavename
        #for i in slave.header.keys():
        #    print i,'\t',slave.header[i]
        
        newslave, footprint = reproject_exact(slave, bheader)
        fits.writeto(outName, newslave, bheader, clobber=True)

        return outName