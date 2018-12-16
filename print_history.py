#!/usr/bin/env python

from astropy.io import fits
import sys

filename=sys.argv[1]

files=fits.open(filename)

heads=files[0].header

if heads['HISTORY']:
	print filename+'\n'
	print heads['HISTORY']
else:
	print '\n*** History is not in header ***\n'