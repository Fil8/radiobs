#!/usr/bin/env python

import sys, os, string
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table
import pyregion
from prettytable import PrettyTable

import montage_wrapper as montage


class makeMask:

	def __init__(self):

		self.rootdir = os.getcwd()+'/'
		self.outmask = self.rootdir+'makeMask_out.fits'

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

	def maskMe(self,fileName,vals,cutoff):
		
		ff=fits.open(fileName)
		dd=ff[0].data
		hh=ff[0].header

		hh = self.cleanHead(hh)
		dd = np.squeeze(dd)


		if len(vals)>0:

			xmin=int(vals[0])
			xmax=int(vals[1])
			ymin=int(vals[2])
			ymax=int(vals[3])
		else:
			xmin=0
			xmax=dd.shape(1)
			ymin=0
			ymax=dd.shape(0)


		new_dd = dd[ymin:ymax,xmin:xmax]
		index = new_dd < cutoff

		new_dd[index==True] = 0.0
		new_dd[index==False] = 1.0

		hh['CRPIX1'] = new_dd.shape[1]/2
		hh['CRPIX2'] = new_dd.shape[0]/2

  
		fits.writeto(self.outmask,new_dd,hh,overwrite=True)

		return 0
#-------------------------------#
#             MAIN              #
#-------------------------------#
masking=makeMask()

# Print help message and exit
if 'help' in sys.argv or '-h' in sys.argv:
	print '\tRun as follows:\n'
	print 'makeMask.py -f <filename.fits> -r <xmin,xmax,ymin,ymax> -c <cutoff>'
	sys.exit()
else: 
	arg=sys.argv

if '-f' in arg:
	fileName=arg[arg.index('-f')+1]
if '-r' in arg: 
	vals=arg[arg.index('-r')+1]
	vals = vals.split(',')
	print vals
if '-c' in arg:
	cutoff = float(arg[arg.index('-c')+1])
else:
	print '\t Parameters missing, run makeMask.py help'
	sys.exit()


masking.maskMe(fileName,vals,cutoff)
	

