#!/usr/bin/env python

from astropy.io import fits
import numpy as np
import sys

filename=sys.argv[1]
outname=sys.argv[2]

print filename
files=fits.open(filename)

datas=files[0].data
heads=files[0].header

datas=np.squeeze(datas)

heads['NAXIS']=2

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

fits.writeto(outname,datas,heads,overwrite=True)