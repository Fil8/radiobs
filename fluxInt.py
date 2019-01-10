#!/usr/bin/env python

__author__ = "Filippo Maccagni"
__copyright__ = "Fil8"
__email__ = "filippo.maccagni@gmail.com"

import sys, os, string
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table
import pyregion
from prettytable import PrettyTable

import montage_wrapper as montage


class flInt:
    

    def __init__(self):

        self.rootdir = os.getcwd()+'/'
        #self.filename = sys.argv[1]

    def openFile(self,filename):
        
        files=fits.open(filename)

        datas = files[0].data
        
        if len(datas.shape)>2:
            datas=np.squeeze(datas)
        
        heads = files[0].header    

        return datas, heads
    
    def readFreq(self,heads,arg):
        
        if '-fr' in arg:
            freq = float(arg[arg.index('-fr')+1])*1e6
        elif 'CRVAL3' in heads:
            freq =  float(heads['CRVAL3'])
        elif 'FREQ' in heads:
            freq = float(heads['FREQ'])

        
        return freq

    def maskData(self,slavename,basename):
        '''
        mask datas within ds9 region estimate noise and background outside of it
        INPUT:
            datas: array of data
            heads: header of file
            mask: array of mask
        OUTPUT:
            masked_data
            background
            noise
            masked_number_of_pixels
        '''
        
        base = fits.open(basename)[0]
        slave = fits.open(slavename)[0]

        headername = self.rootdir+'tmp.hdr'
        outname = string.split(slavename,'.fits')[0]+'_regr.fits'

        print '\t ... reproject image to mask ...'
        if not (slave.header['NAXIS1'] == slave.header['NAXIS1'] and base.header['NAXIS2'] == slave.header['NAXIS2']):
            base.header['BMAJ'] =  slave.header['BMAJ']       
            base.header['BMIN'] =  slave.header['BMIN']
            if slave.header['FREQ']: 
                base.header['FREQ'] =  slave.header['FREQ']    
            elif slave.header['CRVAL3']: 
                base.header['FREQ'] =  slave.header['CRVAL3']    

            fits.writeto(basename,base.data,base.header,overwrite=True)
            montage.mGetHdr(basename,headername)
            montage.mProject(slavename,outname,headername)
            os.remove(headername)
            slave_regr = fits.open(outname)[0]
        else:
            slave_regr = fits.open(slavename)[0]
        
        index = base.data > 0

        #mean stats
        background = np.nanmean(slave_regr.data[index==False])
        self.noise = np.nanstd(slave_regr.data[index==False])
        
        pixels = np.count_nonzero(base.data[index==True])
        noise = np.multiply(self.noise,np.sqrt(pixels))
        
        datas = slave_regr.data[index==True]
        

        return datas, background, noise, pixels, slave_regr.header
 

    def maskDatReg(self,datas,heads,region,cutoff):
        '''
        mask datas within ds9 region estimate noise and background outside of it
        INPUT:
            datas: array of data
            heads: header of file
            region: ds9 region
        OUTPUT:
            masked_data
            background
            noise
            masked_number_of_pixels
        '''

        # set polygonal mask from ds9 region
        r = pyregion.open(region).as_imagecoord(heads)
        shape = (heads['NAXIS2'], heads['NAXIS1'])
        m = r.get_mask(shape=shape)

        #mean stats
        background = np.nanmean(datas[m==False])
        noise = np.nanstd(datas[m==False])
        if cutoff == 0.0:
            self.cutoff = noise*3.
        else:
            self.cutoff = cutoff
            noise=cutoff    

        print self.cutoff

        mm = datas.copy()
        mm[:,:] = 1.
        self.pixels = np.count_nonzero(mm[m==True])
        
        datas[m==False] = np.nan
        index_cut = datas < self.cutoff
        
        datas[index_cut] = np.nan
        

        noise = np.multiply(noise,np.sqrt(self.pixels))
        
        return datas, background, noise, self.pixels

    def cleanHead(self,heads):


        if float(heads['NAXIS']) >2:

            if 'NAXIS3' in heads:
                del heads['NAXIS3']
                del heads['CRVAL3']
                del heads['CDELT3']
                del heads['CRPIX3']
                del heads['CTYPE3']  
                if 'CROTA3' in heads:
                    del heads['CROTA3']
            
        if float(heads['NAXIS']) > 3:
            
            if 'NAXIS4' in heads:
                del heads['NAXIS4']     
                del heads['CRVAL4']
                del heads['CDELT4']
                del heads['CRPIX4']
                del heads['CTYPE4'] 
            if 'CROTA4' in heads:
                del heads['CROTA4']  
        
        heads['NAXIS'] = 2
        return heads

    def measFlux(self,datas,heads,noise,errFlux):

        fluxSum=np.nansum(datas)

        beamArea = 2*np.pi*heads['BMAJ']*3600./2.35482*heads['BMIN']*3600./2.35482
        pixArea = -float(heads['CDELT2']*3600.)*float(heads['CDELT1']*3600.)

        numPixBeam= beamArea/pixArea

        fluxInt = np.divide(fluxSum,numPixBeam)

        return fluxInt,numPixBeam

    def writeTable(self,heads,fluxInt,noise,numPixBeam,freq,errFlux):
        
        noiseInt = np.divide(noise,numPixBeam)
        fluxErr = np.sqrt(np.power(fluxInt/100.*errFlux,2)+np.power(noiseInt,2))          

        alldata = np.array([freq*1e-6,np.round(fluxInt,8),np.round(noiseInt,8),np.round(fluxErr,6),np.round(heads['BMAJ']*3600.,3),
            np.round(heads['BMIN']*3600.,3),np.round(heads['CDELT2']*3600.,0),np.round(numPixBeam,3),np.round(self.cutoff*1e3,4),np.round(self.pixels,0)])

        columnames = ['Frequency [MHz]','Integrated Flux [Jy]','Noise [Jy]', 'Error [Jy]', 'BeamMaj [arcsec]','BeamMin [arcsec]',
                'PixSize [arcsec]','Beam/pix','Flux cutoff [mJy/beam]','PixInt']
       
        self.out_table = self.rootdir+'integratedFluxes.tbl'

        if os.path.exists(self.out_table):
            tt = Table.read(self.out_table, format='ascii')
            tt.add_row(alldata)
        else: 
            tt = Table(alldata,names=columnames,meta={'name': 'Total flux table'})
            
        ascii.write(tt,self.out_table, overwrite=True)

        #print table
        alldata = np.array([freq*1e-6,np.round(fluxInt,5),np.round(noiseInt,6),np.round(fluxErr,6),np.round(heads['BMAJ']*3600.,3),
            np.round(heads['BMIN']*3600.,3),np.round(heads['CDELT2']*3600.,0),np.round(numPixBeam,3),np.round(self.cutoff*1e3,4),np.round(self.pixels,0)])

        t = PrettyTable(columnames)
        t.add_row(alldata)
        
        return t, fluxErr






