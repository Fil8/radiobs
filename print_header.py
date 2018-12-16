#!/usr/bin/env python

from astropy.io import fits
import sys

filename=sys.argv[1]

files=fits.open(filename)

heads=files[0].header

for i in heads.keys():
	print i,'\t',heads[i]

files.close()