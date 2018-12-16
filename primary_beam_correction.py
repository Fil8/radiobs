#!/usr/bin/env python

from astropy.io import fits
from astropy import coordinates as coord
from astropy import wcs
import numpy as np
import os, sys, string



class pbcorr:

	def __init__(self):

		self.rootdir = os.getcwd()+'/'
		self.filename=sys.argv[1]


	def BeamCorrect(self,datas,hdr):

		obs_freq = float(hdr['CRVAL3'])
		

		pb_fwhm = 1.02*(2.99792458E8)/obs_freq/13.5/np.pi*180.
		pb_fwhm_pix = pb_fwhm/hdr['CDELT2']
		x, y = np.meshgrid(np.linspace(-hdr['NAXIS2']/2.,hdr['NAXIS2']/2.,hdr['NAXIS2']), 
						   np.linspace(-hdr['NAXIS1']/2.,hdr['NAXIS1']/2.,hdr['NAXIS1']))
		d = np.sqrt(x*x+y*y)
		sigma, mu = pb_fwhm_pix/2.35482, 0.0
		gauss = np.exp(-( (d-mu)**2 / ( 2.0 * sigma**2 ) ) )
		
		pbcorr = np.divide(datas,gauss)

		outfile=string.split(self.filename,'.fits')[0]
		
		if len(string.split(outfile,'/')) >1:
			outfile=string.split(outfile,'/')[1]
		
		out_pbcorr = self.rootdir + outfile+'_pbcorr.fits'
		out_beam = self.rootdir+'gauss_beam.fits'

		fits.writeto(out_beam,gauss,hdr,overwrite=True)	
		fits.writeto(out_pbcorr,pbcorr,hdr,overwrite=True)


pb = pbcorr()
ff=fits.open(pb.filename)
dats=ff[0].data
heads=ff[0].header
pb.BeamCorrect(dats,heads)


print '''\n\t\t\t------------------------------\n 
			 Primary Beam Correction Done\n
			------------------------------'''