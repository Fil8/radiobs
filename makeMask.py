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
		#self.outmask = self.rootdir+'makeMask_out.fits'

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

	def maskMe(self,fileName,vals,cutoff,options):
		
		ff=fits.open(fileName)
		dd=ff[0].data
		hh=ff[0].header

		hh = self.cleanHead(hh)
		dd = np.squeeze(dd)

		if vals:
			xmin=int(vals[0])
			xmax=int(vals[1])
			ymin=int(vals[2])
			ymax=int(vals[3])
		else:
			xmin=0
			xmax=dd.shape[1]
			ymin=0
			ymax=dd.shape[0]



		new_dd = dd[ymin:ymax,xmin:xmax]
		index = new_dd < cutoff

		new_dd[index==False] = 1.0

		if options == 'nan':
			new_dd[index==True] = np.nan
		elif options == 'zero':
			new_dd[index==True] = 0.0

		hh['CRPIX1'] = new_dd.shape[1]/2
		hh['CRPIX2'] = new_dd.shape[0]/2

		outfile=string.split(fileName,'.fits')[0]
		outfile = outfile+'_mask.fits'

  
		fits.writeto(outfile,new_dd,hh,overwrite=True)
		
				


		print '''\n\t\t\t-----------------------------
			   Mask from cutoff done
			-----------------------------\n'''

		return 0

	def cutInCentre(self,filename,pixSizeX,pixSizeY):

		ff=fits.open(fileName)
		dd=ff[0].data
		hh=ff[0].header

		hh = self.cleanHead(hh)
		dd = np.squeeze(dd)

		xmin = int(np.round(hh['CRPIX1'],0)-np.round(pixSizeX/2.,0))
		xmax = int(np.round(hh['CRPIX1'],0)+np.round(pixSizeX/2.,0))

		ymin = int(np.round(hh['CRPIX2'],0)-np.round(pixSizeY/2.,0))
		ymax = int(np.round(hh['CRPIX2'],0)+np.round(pixSizeY/2.,0))
		
		newDD = dd[ymin:ymax,xmin:xmax]
		naxis1 = newDD.shape[1]
		naxis2 = newDD.shape[0]
		crpix1 = newDD.shape[1]/2
		crpix2 = newDD.shape[0]/2

		hh['CRPIX1'] = crpix1
		hh['CRPIX2'] = crpix2
		hh['NAXIS1'] = naxis1
		hh['NAXIS2'] = naxis2

		outfile=string.split(fileName,'.fits')[0]
		outfile = outfile+'_cutCtr.fits'

		fits.writeto(outfile,newDD,hh,overwrite=True)

		print '''\n\t\t\t-----------------------------
			   Cutout from Centre done
			-----------------------------\n'''

		return 0

	def to32Bits(self,filename):

		ff=fits.open(filename)
		dd=ff[0].data
		hh=ff[0].header
		dd=dd.astype('float32')

		outfile=string.split(fileName,'.fits')[0]
		outfile = outfile+'_bt32.fits'

		fits.writeto(outfile,dd,hh,overwrite=True)

		print '''\n\t\t\t-----------------------------
			   Conversion to 32-BIT done
			-----------------------------\n'''
#-------------------------------#
#             MAIN              #
#-------------------------------#
masking=makeMask()

# Print help message and exit
if 'help' in sys.argv or '-h' in sys.argv:
	print '\tTo make mask run as follows:\n'
	print 'makeMask.py -m -f <filename.fits> -r <xmin,xmax,ymin,ymax> -c <cutoff> -o <nan/zero>'
	print '\tTo make mask run as follows:\n'
	print 'makeMask.py -cutCtr -f <filename.fits> -px X-size_pixel -py Y-size-pixel'
	print '\tTo make mask run as follows:\n'
	print 'makeMask.py -cv32 -f <filename.fits>'
	sys.exit()
else: 
	arg=sys.argv

if '-f' in arg:
	fileName=arg[arg.index('-f')+1]
if '-r' in arg: 
	vals=arg[arg.index('-r')+1]
	vals = vals.split(',')
else: 
	vals = None

if '-o' in arg:
	options = arg[arg.index('-o')+1]
if '-c' in arg:
	cutoff = float(arg[arg.index('-c')+1])
if '-px' in arg:
	pixX = float(arg[arg.index('-px')+1])
if '-py' in arg:
	pixY = float(arg[arg.index('-py')+1])
#else:
#	print '\t Parameters missing, run makeMask.py help'
#	sys.exit()

if '-m' in arg:
	masking.maskMe(fileName,vals,cutoff,options)
elif '-cutCtr' in arg:
	masking.cutInCentre(fileName,pixX,pixY)
elif ('-cv32' in arg) and (fileName!=None):
	masking.to32Bits(fileName)

