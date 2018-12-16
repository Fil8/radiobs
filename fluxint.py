#!/usr/bin/env python

import sys
import numpy as np
from astropy.io import fits
import pyregion
from prettytable import PrettyTable

def MaskDatReg(datas,heads,region,cutoff):
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
    print cutoff
    index_cut = datas < cutoff
    
    datas[index_cut] = np.nan

    noise = np.multiply(noise,np.sqrt(pixels))
    
    return datas, background, noise, pixels



def CleanHead(heads):

    heads['NAXIS'] = 2

    if float(heads['NAXIS']) >2:
        del heads['NAXIS3']
        del heads['CRVAL3']
        del heads['CDELT3']
        del heads['CRPIX3']
        del heads['CTYPE3']  

    if float(heads['NAXIS']) > 3:
        heads['NAXIS'] = 2

        del heads['NAXIS3']
        del heads['CRVAL3']
        if 'CDELT3' in heads:
            del heads['CDELT3']
        del heads['CRPIX3']
        del heads['CTYPE3']  
        del heads['NAXIS4']     
        del heads['CRVAL4']
        if 'CDELT4' in heads:
            del heads['CDELT4']
        del heads['CRPIX4']
        del heads['CTYPE4'] 
        if 'CROTA3' in heads:
            del heads['CROTA3']
        if 'CROTA4' in heads:
            del heads['CROTA4']  

    return heads

def MeasFlux(datas,heads):


    fluxsum=np.nansum(maskedat)

    beam_area = 2*np.pi*heads['BMAJ']*3600./2.35482*heads['BMIN']*3600./2.35482
    pix_area = -float(heads['CDELT2']*3600.)*float(heads['CDELT1']*3600.)

    number_pix_beam= beam_area/pix_area

    fluxint = np.divide(fluxsum,number_pix_beam)


    return fluxint,number_pix_beam

filename=sys.argv[1]
region=sys.argv[2]
tableout=sys.argv[3]
cutoff=float(sys.argv[4])


files=fits.open(filename)

datas = files[0].data
datas=np.squeeze(datas)
heads = files[0].header    

freqstart = heads['CRVAL3']

heads=CleanHead(heads)    

maskedat, background, noise, pixels=MaskDatReg(datas,heads,region,cutoff)

fluxint,number_pix_beam=MeasFlux(maskedat,heads)
noiseint = np.divide(noise,number_pix_beam)

alldata = np.array([freqstart*1e-6,np.round(fluxint,3),np.round(noiseint,5),np.round(heads['BMAJ']*3600.,3),np.round(heads['BMIN']*3600.,3),np.round(heads['CDELT2']*3600.,0),np.round(number_pix_beam,3),np.round(pixels,0)])
columnames = ['Frequency [MHz]','Integrated Flux [Jy]','Noise [Jy]','BeamMaj [arcsec]','BeamMin [arcsec]', 'PixSize [arcsec]',
              'Beam/pix','PixInt']
t = PrettyTable(columnames)
t.add_row(alldata)

print t