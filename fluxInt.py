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

#import montage_wrapper as montage


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

        heads = self.cleanHead(heads)

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
 

    def maskDatReg(self,dats,heads,region,cutoff):
        '''
        mask datas within ds9 region estimate noise and background outside of it
        INPUT:
            datas: array of data
            heads: header of file
            region: ds9 region
        OUTPUT:
            masked_data
            background
            rms
            masked_number_of_pixels
        '''
        datas = dats.copy()
        # set polygonal mask from ds9 region
        r = pyregion.open(region).as_imagecoord(heads)
        shape = (heads['NAXIS2'], heads['NAXIS1'])
        m = r.get_mask(shape=shape)

        #mean stats
        background = np.nanmean(datas[m==False])
        noise = np.nanstd(datas[m==False])
        
        #if cutoff == 0.0:
        self.cutoff = np.nan
        #elif cutoff < 0:
        #    self.cutoff = noise
        #else:
        #    self.cutoff = cutoff
        #    noise=cutoff            

        datas[m==False] = np.nan
        index_cut = datas < self.cutoff
        
        datas[index_cut] = np.nan
        
        mm = datas.copy()
        mm[:,:] = 1.
        self.pixels = np.count_nonzero(mm[m==True])
        
        #noise = np.multiply(noise,np.sqrt(self.pixels))
        
        return datas, background, noise, self.pixels

    def noiseReg(self,datas,heads,region):
        '''
        mask datas within ds9 region estimate noise and background outside of it
        INPUT:
            datas: array of data
            heads: header of file
            region: ds9 region
        OUTPUT:
            background
            rms
            number_of_pixels_in_region
        '''

        # set polygonal mask from ds9 region
        r = pyregion.open(region).as_imagecoord(heads)
        shape = (heads['NAXIS2'], heads['NAXIS1'])
        m = r.get_mask(shape=shape)

        #mean stats
        background = np.nanmean(datas[m==False])
        noise = np.nanstd(datas[m==False]) 

        mm = datas.copy()
        mm[np.isnan(mm)] = 0.
        pixels = np.count_nonzero(mm[m==True])
        return background, noise, pixels

    def noiseMultiReg(self,ldata,lhead,region_dir):
        
        r = [f for f in os.listdir(region_dir) if os.path.isfile(os.path.join(region_dir, f))]

        mean_values = np.zeros(len(r))
        med_values = np.zeros(len(r))

        flux_values = np.zeros(len(r))
        pixels = 0.
        for i in xrange(0,len(r)):
            region_name = region_dir+r[i]
            #print region_name
            re = pyregion.open(region_name).as_imagecoord(lhead)
            shape = (lhead['NAXIS2'], lhead['NAXIS1'])
            m_noise = re.get_mask(shape=shape)
            mask_tmp = np.copy(ldata)
            #print np.nansum(mask_tmp[m_noise==True])

            #mask_tmp[m==True] = np.nan
            flux_values[i] = np.nansum(mask_tmp[m_noise==True])
            mean_values[i] = np.nanmean(mask_tmp[m_noise==True])
            med_values[i] = np.nanmedian(mask_tmp[m_noise==True])

            pixels += np.count_nonzero(mask_tmp[m_noise==True])
        
        noise = np.nanstd(flux_values)
        back = np.nanmean(mean_values)
        backmed = np.nanmedian(med_values)
        perc = np.percentile(flux_values,67)

        pixels = pixels/len(r)

        return back,noise, pixels

    def cleanHead(self,heads):
        #base = fits.open(fileName)

        #heads = base[0].header
        #datas = base[0].data
        if 'NAXIS3' in heads:
            del heads['NAXIS3']
        if 'CRVAL3' in heads:
            del heads['CRVAL3']
        if 'CDELT3' in heads:
            del heads['CDELT3']
        if 'CRPIX3' in heads: 
            del heads['CRPIX3']
        if 'CTYPE3' in heads:        
            del heads['CTYPE3']  
        if 'CROTA3' in heads:
            del heads['CROTA3']
        if 'CUNIT3' in heads:
            del heads['CUNIT3']        
            
        if 'NAXIS4' in heads:
            del heads['NAXIS4']     
        if 'CRVAL4' in heads:        
            del heads['CRVAL4']
        if 'CDELT4' in heads:    
            del heads['CDELT4']
        if 'CRPIX4' in heads:    
            del heads['CRPIX4']
        if 'CTYPE4' in heads:    
            del heads['CTYPE4'] 
        if 'CROTA4' in heads:
            del heads['CROTA4']  
        if 'CUNIT4' in heads:
            del heads['CUNIT4']
        
        if 'WCSAXES' in heads:
            heads['WCSAXES'] = 2
        
        if 'PC1_1' in heads:   
            del heads['PC1_1']            
        if 'PC2_1' in heads:   
            del heads['PC2_1']           
        if 'PC3_1' in heads:   
            del heads['PC3_1']
        if 'PC4_1' in heads:   
            del heads['PC4_1']    
        if 'PC1_2' in heads:
            del heads['PC1_2']
        if 'PC2_2' in heads:
            del heads['PC2_2']
        if 'PC3_2' in heads:
            del heads['PC3_2']
        if 'PC4_2' in heads:
            del heads['PC4_2']
        if 'PC1_3' in heads:
            del heads['PC1_3']
        if 'PC2_3' in heads:
            del heads['PC2_3']
        if 'PC3_3' in heads:    
            del heads['PC3_3']
        if 'PC4_3' in heads:    
            del heads['PC4_3']
        if 'PC1_4' in heads:
            del heads['PC1_4']            
        if 'PC2_4' in heads:
            del heads['PC2_4']
        if 'PC3_4' in heads:
            del heads['PC3_4']
        if 'PC4_4' in heads:
            del heads['PC4_4']
        if 'PV2_2' in heads:
            del heads['PV2_1']
        if 'PV2_2' in heads:
            del heads['PV2_2']

        heads['NAXIS'] = 2

        #fits.writeto(fileName,datas,heads,overwrite=True)

        return heads

    def measFlux(self,datas,heads,errFlux,option):

        fluxSum=np.nansum(datas)
        if option=='RgCv':
            heads['BMAJ'] = 18.5/3600.
            heads['BMIN'] = 9./3600.

        beamArea = 2*np.pi*float(heads['BMAJ'])*3600./2.35482*float(heads['BMIN'])*3600./2.35482
        print heads['BMAJ']
        pixArea = -float(heads['CDELT2']*3600.)*float(heads['CDELT1']*3600.)

        numPixBeam= beamArea/pixArea


        fluxInt = np.divide(fluxSum,numPixBeam)

        return fluxInt,numPixBeam

    def writeTable(self,heads,fluxInt,noise,numPixBeam,freq,errFlux,outTable='integratedFluxes.tbl'):
        
        #noiseInt = np.divide(noise,numPixBeam)
        noiseInt = float(noise)
        fluxErr = np.sqrt(np.power(fluxInt/100.*errFlux,2)+np.power(noiseInt,2))          
        alldata = np.array([freq*1e-6,np.round(fluxInt,8),np.round(noiseInt,8),np.round(fluxErr,6),np.round(float(heads['BMAJ'])*3600.,3),
            np.round(float(heads['BMIN'])*3600.,3),np.round(float(heads['CDELT2'])*3600.,0),np.round(numPixBeam,3),np.round(self.cutoff*1e3,4),np.round(self.pixels,0)])

        columnames = ['Frequency [MHz]','Integrated Flux [Jy]','Noise [Jy]', 'Error [Jy]', 'BeamMaj [arcsec]','BeamMin [arcsec]',
                'PixSize [arcsec]','Beam/pix','Flux cutoff [mJy/beam]','PixInt']
       
        #self.out_table = self.rootdir+outTable

        if os.path.exists(outTable):
            tt = Table.read(outTable, format='ascii')
            tt.add_row(alldata)
        else: 
            tt = Table(alldata,names=columnames,meta={'name': 'Total flux table'})
            
        ascii.write(tt,outTable, overwrite=True)

        #print table
        alldata = np.array([freq*1e-6,np.round(fluxInt,5),np.round(noiseInt,6),np.round(fluxErr,6),np.round(float(heads['BMAJ'])*3600.,3),
            np.round(float(heads['BMIN'])*3600.,3),np.round(float(heads['CDELT2'])*3600.,0),np.round(numPixBeam,3),np.round(self.cutoff*1e3,4),np.round(self.pixels,0)])

        t = PrettyTable(columnames)
        t.add_row(alldata)
        
        return t, fluxErr










