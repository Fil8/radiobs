#!/usr/bin/env python

import sys, os
from astropy.io import fits
import numpy as np


class chHeader:


	def __init__(self):

		self.rootdir = os.getcwd()+'/'
		#self.filename = sys.argv[1]

	def wscleanChHead(self,fileName,outName):

		files=fits.open(fileName)

		datas=files[0].data
		heads=files[0].header

		datas=np.squeeze(datas)
		

		heads['NAXIS']=2

		try :
			heads['FREQ']
		except KeyError:
			heads['FREQ'] = heads['CRVAL3']
		else:
			pass
		if heads['NAXIS3']:
			del heads['NAXIS3']

		if heads['NAXIS4']:
			del heads['NAXIS4']

		if heads['CRVAL4']:
			del heads['CRVAL4']

		if heads['CRVAL3']:
			del heads['CRVAL3']

		if heads['CRPIX4']:
			del heads['CRPIX4']

		if heads['CRPIX3']:
			del heads['CRPIX3']

		if heads['CDELT3']:
			del heads['CDELT3']

		if heads['CDELT4']:
			del heads['CDELT4']

		if heads['CTYPE4']:
			del heads['CTYPE4']

		if heads['CTYPE3']:
			del heads['CTYPE3']

		fits.writeto(outName,datas,heads,overwrite=True)

		return 0

	def putHead(self,fileName,key,value):

		files=fits.open(fileName)

		datas=files[0].data
		heads=files[0].header

		heads[key] = value

		fits.writeto(fileName,datas,heads,overwrite=True)

		return 0 

#-------------------------------#
#             MAIN              #
#-------------------------------#
chHead=chHeader()

# Print help message and exit
if 'help' in sys.argv or '-h' in sys.argv:
	print 'Run as follows'
	print 'chHead.py -f function: <wsclHead OR putHead> -i <filename.fits> OPTIONS -o <filename.fits> -key <header_key> -val <header_value>'
	sys.exit()
else: arg=sys.argv

if '-i' in arg:
	fileName=arg[arg.index('-i')+1]
else:
	print 'Run as follows'
	print 'chHead.py -f function: <wsclHead OR putHead> -i <filename.fits> OPTIONS -o <filename.fits>'

if '-o' in arg:
	outName=arg[arg.index('-o')+1]
else:
	outName = fileName

if '-f' in arg:
	function = arg[arg.index('-f')+1]
else:
	print 'Run as follows'
	print 'chHead.py -f function: <wsclHead OR putHead> -i <filename.fits> OPTIONS -o <filename.fits>'

if function == 'wsclHead':
	chHead.wscleanChHead(fileName,outName)


if function == 'putHead':

	key = arg[arg.index('-key')+1]
	value = arg[arg.index('-val')+1]
	chHead.putHead(fileName,key,value)

print '''\t+---------+\n\t Head done\n\t+---------+'''
