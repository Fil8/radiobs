#!/usr/bin/env python

__author__ = "Filippo Maccagni"
__copyright__ = "Fil8"
__email__ = "filippo.maccagni@gmail.com"

import sys, os, string
import numpy as np
from astropy.io import fits, ascii
from scipy import optimize

class fitSrc:
    

    def __init__(self):

        self.rootdir = os.getcwd()+'/'
        #self.filename = sys.argv[1]
    
    def gaussian(self,height, center_x, center_y, width_x, width_y):
        """Returns a gaussian function with the given parameters"""
        width_x = float(width_x)
        width_y = float(width_y)
        return lambda x,y: height*np.exp(
                    -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

    def moments(self,data):
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

    def fitgaussian(self,data):
        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution found by a fit"""
        params = self.moments(data)
        errorfunction = lambda p: np.ravel(self.gaussian(*p)(*np.indices(data.shape)) -
                                     data)
        p, success = optimize.leastsq(errorfunction, params)
        return p


    def gaus2Dfit(self,data):
        
        data = np.squeeze(data)
        pars = self.fitgaussian(data)

        return pars