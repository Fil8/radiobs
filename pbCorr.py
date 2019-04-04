#!/usr/bin/env python

from astropy.io import fits
from astropy import coordinates as coord
from astropy import wcs
import numpy as np
import os, sys, string



class pbcorr:

	def __init__(self):

		self.rootdir = os.getcwd()+'/'


	def BeamCorrect(self,fileName,datas,hdr,telescope):

		try :
			hdr['CRVAL3']
		except KeyError:
			obs_freq = float(hdr['FREQ'])	
		else:
			obs_freq = float(hdr['CRVAL3'])

		if 'RESTFREQ' in hdr:
			hdr['FREQ'] = hdr['RESTFREQ']

		if telescope == 'MeerKAT':
			ant = 13.5
		elif telescope == 'VLA':
			ant = 25.
		elif telescope == 'ACA':
			ant = 7.

		pb_fwhm = 1.02*(2.99792458E8)/obs_freq/ant/np.pi*180.
		pb_fwhm_pix = pb_fwhm/hdr['CDELT2']
		x, y = np.meshgrid(np.linspace(-hdr['NAXIS1']/2.,hdr['NAXIS1']/2.,hdr['NAXIS1']), 
						   np.linspace(-hdr['NAXIS2']/2.,hdr['NAXIS2']/2.,hdr['NAXIS2']))
		d = np.sqrt(x*x+y*y)
		sigma, mu = pb_fwhm_pix/2.35482, 0.0
		gauss = np.exp(-( (d-mu)**2 / ( 2.0 * sigma**2 ) ) )
		
		pbcorr = np.divide(datas,gauss)
		outfile=string.split(fileName,'.fits')[0]
		#if len(string.split(outfile,'/')) >1:
	#		outfile=string.split(outfile,'/')[1]
		
		out_pbcorr = outfile+'_pbcorr.fits'
		out_beam = outfile+'_gauss_beam.fits'

		fits.writeto(out_beam,gauss,hdr,overwrite=True)	
		fits.writeto(out_pbcorr,pbcorr,hdr,overwrite=True)


		return out_pbcorr
