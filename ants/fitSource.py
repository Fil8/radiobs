#!/usr/bin/env python

import sys, os, string
import numpy as np
from astropy.io import fits, ascii
from scipy import optimize
from scipy import fftpack

from astropy.table import Table
from prettytable import PrettyTable

class fitsrc:
    

    def __init__(self):

        self.rootdir = os.getcwd()+'/'

    def gauss1D(self,x, *p):
        
        a, b, c, d = p
        y = a*np.exp(-np.power((x - b), 2.)/(2. * c**2.)) + d

        return y

    def convolve2DGauss(self, star, psf):
        star_fft = fftpack.fftshift(fftpack.fftn(star))
        psf_fft = fftpack.fftshift(fftpack.fftn(psf))
        return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(star_fft*psf_fft)))

    def deconvolve2DGauss(self, star, psf):
        star_fft = fftpack.fftshift(fftpack.fftn(star))
        psf_fft = fftpack.fftshift(fftpack.fftn(psf))
        return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(star_fft/psf_fft)))
    
    def gaussian2D(self,height, center_x, center_y, width_x, width_y):
        """Returns a gaussian function with the given parameters"""
        width_x = float(width_x)
        width_y = float(width_y)
        return lambda x,y: height*np.exp(
                    -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

    def moments2D(self,data):
        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution by calculating its
        moments """
        total = data.sum()
        X, Y = np.indices(data.shape)
        x = (X*data).sum()/total
        y = (Y*data).sum()/total
        col = data[:, int(y)]
        width_x = np.sqrt(np.abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
        row = data[int(x), :]
        width_y = np.sqrt(np.abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
        height = data.max()

        return height, x, y, width_x, width_y


    def fit2Dgaussian(self,data):
        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution found by a fit"""
        params = self.moments2D(data)
        errorfunction = lambda p: np.ravel(self.gaussian2D(*p)(*np.indices(data.shape)) -
                                     data)
        p, success = optimize.leastsq(errorfunction, params)
        return p

    def fit1Dgaussian(self,data):
        """Returns (height, x, width_x)
        the gaussian parameters of a 1D distribution found by a fit"""
        
        X = np.arange(data.size)
        x = np.sum(X*data)/np.sum(data)
        
        width = np.sqrt(np.abs(np.sum((X-x)**2*data)/np.sum(data)))

        max = data.max()

        fit = lambda t : max*np.exp(-(t-x)**2/(2*width**2))
        
        return max,x,width

    def gaus2Dfit(self,data):
        
        data = np.squeeze(data)
        pars = self.fit2Dgaussian(data)
        return pars


    def gaus1Dfit(self,data):

        data = np.squeeze(data)
        pars = self.fit1Dgaussian(data)
        return pars

    def writeGaus2DfitTable(self,region,outTable,pars):
        
        pars = pars.astype('str')
        alldata = np.insert(pars,0,region)

        columnames = ['RegionName','IntFlux','PeakIm','PeakGaus','centre X','centre Y', 'width X', 'width Y']
       
        tt = Table(alldata,names=columnames,meta={'name': 'Jet Counter-Jet table'})
       
        if os.path.exists(outTable):
            tt = Table.read(outTable, format='ascii')
            tt.add_row(alldata)
        else: 
            tt = Table(alldata,names=columnames,meta={'name': 'Total flux table'})    
        ascii.write(tt,outTable, overwrite=True)
        
        #print alldata
 
        t = PrettyTable(columnames)
        t.add_row(alldata)
        
    def writeGaus1DfitTable(self,region,outTable,pars):
        
        pars = pars.astype('str')
        alldata = np.insert(pars,0,region)

        columnames = ['RegionName','IntFlux','PeakIm','PeakGaus','centre X', 'width X']
       
        tt = Table(alldata,names=columnames,meta={'name': 'Jet Counter-Jet table'})
       
        if os.path.exists(outTable):
            tt = Table.read(outTable, format='ascii')
            tt.add_row(alldata)
        else: 
            tt = Table(alldata,names=columnames,meta={'name': 'Total flux table'})    
        ascii.write(tt,outTable, overwrite=True)
        
        #print alldata
 
        t = PrettyTable(columnames)
        t.add_row(alldata)
            

        return t       