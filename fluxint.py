#!/usr/bin/env python

import sys, os, string
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table
import pyregion
from prettytable import PrettyTable

import montage_wrapper as montage


class fluxint:
    

    def __init__(self):

        self.rootdir = os.getcwd()+'/'
        #self.filename = sys.argv[1]
        self.out_table = self. rootdir+'integrated_flux.tbl'

    def openFile(self,filename):
        
        files=fits.open(filename)

        datas = files[0].data
        
        if len(datas.shape)>2:
            datas=np.squeeze(datas)
        
        heads = files[0].header    

        return datas, heads
    
    def readFreq(self,heads,arg):
        
        if 'CRVAL3' in heads:
            freq = heads['CRVAL3']
        elif '-fr' in arg:
            freq = float(arg[arg.index('-fr')+1])*1e6

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
        
            montage.mGetHdr(basename,headername)
            montage.mProject(slavename,outname,headername)
            os.remove(headername)
            slave_regr = fits.open(outname)[0]
        else:
            slave_regr = fits.open(slavename)[0]

        print '\t ... mask image ...'
        
        index = base.data > 0

        #mean stats
        background = np.nanmean(slave_regr.data[index==False])
        noise = np.nanstd(slave_regr.data[index==False])
        
        pixels = np.count_nonzero(base.data[index==True])
        
        noise = np.multiply(noise,np.sqrt(pixels))
        
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
        
        mm = datas.copy()
        mm[:,:] = 1.
        pixels = np.count_nonzero(mm[m==True])
        
        datas[m==False] = np.nan
        index_cut = datas < cutoff
        
        datas[index_cut] = np.nan
        

        noise = np.multiply(noise,np.sqrt(pixels))
        
        return datas, background, noise, pixels



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

    def measFlux(self,datas,heads):

        print '\t ... measure flux ...'

        fluxsum=np.nansum(datas)

        beam_area = 2*np.pi*heads['BMAJ']*3600./2.35482*heads['BMIN']*3600./2.35482
        pix_area = -float(heads['CDELT2']*3600.)*float(heads['CDELT1']*3600.)

        number_pix_beam= beam_area/pix_area

        fluxint = np.divide(fluxsum,number_pix_beam)

        return fluxint,number_pix_beam

    def writeTable(self,heads,fluxint,noise_pix_beam,freq):
        
        print '\t ... write table ...'

        noiseint = np.divide(noise,number_pix_beam)
        
        alldata = np.array([freq*1e-6,np.round(fluxint,8),np.round(noiseint,8),np.round(heads['BMAJ']*3600.,3),
            np.round(heads['BMIN']*3600.,3),np.round(heads['CDELT2']*3600.,0),np.round(number_pix_beam,3),np.round(pixels,0)])

        columnames = ['Frequency [MHz]','Integrated Flux [Jy]','Noise [Jy]','BeamMaj [arcsec]','BeamMin [arcsec]',
                'PixSize [arcsec]','Beam/pix','PixInt']
        
        if os.path.exists(self.out_table):
            tt = Table.read(self.out_table, format='ascii')
            tt.add_row(alldata)
        else: 
            tt = Table(alldata,names=columnames,meta={'name': 'Total flux table'})
            
        ascii.write(tt,self.out_table, overwrite=True)

        #print table
        alldata = np.array([freq*1e-6,np.round(fluxint,3),np.round(noiseint,5),np.round(heads['BMAJ']*3600.,3),
            np.round(heads['BMIN']*3600.,3),np.round(heads['CDELT2']*3600.,0),np.round(number_pix_beam,3),np.round(pixels,0)])

        t = PrettyTable(columnames)
        t.add_row(alldata)
        print t
        
        return 0




#-------------------------------#
#             MAIN              #
#-------------------------------#
fint=fluxint()

# Print help message and exit
if 'help' in sys.argv or '-h' in sys.argv:
    print 'Run as follows'
    print 'fluxint.py -f <filename.fits> OR -t <filename.tbl> AND -m <maskname.fits> OR -r <region.fits> OPTIONS (-c <cutoff> -fr <observed_frequency in MHz>)'
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
    print region

if '-c' in arg:
    cutoff =float(arg[arg.index('-c')+1])
else:
    cutoff = 1e-3


if ('-f' in arg and '-r' in arg):
    print 'ciao'
    print '\t Executing File+Region Combo' 
    datas,heads = fint.openFile(fileName)
    datas=np.squeeze(datas)
    heads = heads=fint.cleanHead(heads)
    maskedData, background, noise, pixels=fint.maskDatReg(datas,heads,region,cutoff)

    freq=fint.readFreq(heads,arg)
    fluxint, number_pix_beam=fint.measFlux(maskedData,heads)
    fint.writeTable(heads,fluxint,number_pix_beam,freq)

if ('-f' in arg and '-m' in arg):
    print '\t Executing File+Mask Combo' 
    maskedData, background, noise, pixels, newHeads = fint.maskData(fileName,maskName)
    heads=newHeads

    freq=fint.readFreq(heads,arg)
    fluxint, number_pix_beam=fint.measFlux(maskedData,heads)
    fint.writeTable(heads,fluxint,number_pix_beam,freq)

if ('-t' in arg and '-m' in arg):
    print '\t Executing FileList+Mask Combo'
    for i in xrange(0,len(tableFileNames.columns[0])):

        fileName = tableFileNames.columns[0][i]
        maskedData, background, noise, pixels, newHeads = fint.maskData(fileName,maskName)
        heads=newHeads
        
        if 'Frequency' in tableFileNames.dtype.names:
            freq = tableFileNames.columns[1][i]*1e6

        fluxint, number_pix_beam=fint.measFlux(maskedData,heads)

        fint.writeTable(heads,fluxint,number_pix_beam,freq)
