#!/usr/bin/env python

__author__ = "Filippo Maccagni"
__copyright__ = "Fil8"
__email__ = "filippo.maccagni@gmail.com"

import sys, os, string
import numpy as np
from astropy.io import fits, ascii
from scipy import optimize
#from scipy.optimize import curve_fit
from scipy import fftpack

from astropy.table import Table

from prettytable import PrettyTable

class fitSrc:
    

    def __init__(self):

        self.rootdir = os.getcwd()+'/'
        #self.filename = sys.argv[1]
    
    def gauss1D(self, x, amp, xc, sigma):
    
        return amp*np.exp( -(x-xc)**2 / (2*sigma**2)) / np.sqrt(2*np.pi*sigma**2)
    
    def gauss2D(self,(x,y), amplitude, x0, y0, sigma_x, sigma_y, theta, offset):
        
        """Returns a gaussian function with the given parameters"""
        #width_x = float(width_x)
        #width_y = float(width_y)

        #g= amplitude*np.exp(
        #            -(((x0-x)/sigma_x)**2+((y0-y)/sigma_y)**2)/2)        
        
        #(x, y) = xdata_tuple                                                        
        xo = float(x0)                                                              
        yo = float(y0)                                                              
        a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)   
        b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)    
        c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)   
        g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)         
                        + c*((y-yo)**2)))                                   
        
        return g.ravel()   
    
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

    def fitgauss1D(self,data):
        p_initial = [np.nanmax(data), data.shape[0]/2., 3.]
        

        # Create our data sets. Perturb the y-data with randomness and
        # generate completely random data for the errors.
        x = np.arange(data.shape[0])
        # Use curve_fit to fit the gauss function to our data. Use the
        # unperturbed p_initial as our initial guess.
        popt, pcov = optimize.curve_fit(self.gauss1D, x, data, p_initial)
        errors =  np.sqrt(np.diag(pcov))
        return popt

    def fitgauss2D(self,data):
        print data.shape
        
        ind = np.unravel_index(np.argmax(data, axis=None), data.shape)
        print ind
        print np.nanmax(data)
        p_initial = [np.nanmax(data), ind[0], ind[1], 3., 3., 10., 0.]
        
        # Create our data sets. Perturb the y-data with randomness and
        # generate completely random data for the errors.
        y = np.arange(data.shape[0])
        x = np.arange(data.shape[1])
        x, y = np.meshgrid(x, y)

        xdata = np.vstack((x.ravel(),y.ravel()))
        data = data.ravel()
        # Use curve_fit to fit the gauss function to our data. Use the
        # unperturbed p_initial as our initial guess.
        popt, pcov = optimize.curve_fit(self.gauss2D, xdata, data, p_initial)
        errors =  np.sqrt(np.diag(pcov))
        return popt


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