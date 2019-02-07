#!/usr/bin/env python

from astropy.io import fits
import pbCorr
import numpy as np
import os, sys, string

pb = pbCorr.pbcorr()

fileName=sys.argv[1]
telescope=sys.argv[2]
ff=fits.open(fileName)
dats=ff[0].data
heads=ff[0].header
dats=np.squeeze(dats)
pb.BeamCorrect(fileName,dats,heads,telescope)


print '''\n\t\t\t------------------------------\n 
			 Primary Beam Correction Done\n
			------------------------------'''
