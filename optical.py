__author__ = "Filippo Maccagni"
__copyright__ = "Fil8"
__email__ = "filippo.maccagni@gmail.com"

import sys,string,os,math
import numpy as np
from radiobs import Cosmo
from astropy.io import fits
from photutils.isophote import EllipseGeometry, Ellipse, build_ellipse_model
from photutils.isophote.sample import CentralEllipseSample
from photutils.isophote.fitter import CentralEllipseFitter
from photutils import EllipticalAperture


class optical:
	
	def __init__(self):
		#define constants
		RAD2DEG=180./math.pi
		HI=1.42040575177e+09 #Hz
		TSPIN=100            #K
		MSUN=1.98855e33      #g
		MHI=1.6749E-24       #g
		CHI=2.36E5
		PC=3.08567758E18    #cm
		JANSKY=1e-23        #erg/scm2Hz
		C=2.99792458E10     #cm/s
		G=6.6742E-08        #cm3kg-1s-1      
		MP=1.67492728E-24   #g
		SIGMAT=6.66524E-25  #cm2


	def abs_mag(self,m,z):
		'''
		Estimates the absolute magnitude of a source give its redshift and apparent magnitude
			INPUT
				z: redshift of the source (z)
				m: apparent magnitude of the source
			OUTPUT
				M: absolute magnitude of the source
		'''			

		dl=c.lum_dist(z)
		M=m-5*(math.log10(dl/3.085678e18)-1)
		
		return M

	def vflat_tf(self,M,z):
		'''
		Estimates the flat rotational velocity (large radii) of the source according to the 
		Tully-Fisher relation
		INPUT
				M: absolute magnitude of the source
					output of observe.abs_mag
			OUTPUT
				M: absolute magnitude of the source
		'''

		vflattf= 10**((M-self.m26)/self.s_tully+2.6)/2	#km/s

		print 'Rotational velocity (TullyFisher flat) = '+str(round(vflattf,3))+' km/s'
	
		return vflattf

	def iso_fit(self,infits,outdir,geometry,Norm):
		'''
		Isophotal fitting of a .fits photometric image
			INPUT
				infits: photometric image to fit in .fits format
				outdir: directory where to save fit results
				geometry: geometry to fit (EllipseGeometry module)
				Norm: True/False - normalize fluxes of image between 0 and 1
			OUTPUT
				model.fits: image of the fitted model
				residual.fits: image of the residuals between the model and the data
				iso_table.txt: isolist results converted to table
		'''	

		#create output directory
		if os.path.isdir(outdir) == False:
			os.mkdir(outdir)

		# open input file
		inlist = fits.open(infits)
		self.datain = inlist[0].data
		#normalize data between 0 and 1
		if Norm == True:
			self.data = (self.datain-np.min(self.datain))/(np.max(self.datain)-np.min(self.datain))
		else:
			self.data = self.datain
		self.head = inlist[0].header

		#create instance of Ellipse, inputting the data to be fitted and the initial ellipse geometry object
		ellipse = Ellipse(self.data, geometry)
		# fit data
		#isolist = ellipse.fit_image(maxsma=(self.head['NAXIS1']/2.8))
		isolist = ellipse.fit_image(maxsma=(self.head['NAXIS1']/2.8))

		# -------------------------- #
		# Normalize the central peak
		# we use the find_center method to fine tune and confirm our
		# choice of central pixel. 
		#geometry.find_center(self.data)

		# then we build a CentralEllipseSample instance and call the 
		# CentralEllipseFitter on it.
		#sample = CentralEllipseSample(self.data, 0., geometry=geometry)
		#fitter = CentralEllipseFitter(sample)
		#center = fitter.fit()
		#del isolist[0]
		#isolist.append(center)
		#isolist.sort()
		# -------------------------- #


		#isolist.fix_geometry(isolist.get_closest(isolist.sma[-1]))


		#tab = isolist.to_table()

		# build an elliptical model image from the IsophoteList object 
		model_image = build_ellipse_model(self.data.shape, isolist)

		#if Norm == True:
		#	model_image = model_image*(np.max(self.datain)-np.min(self.datain)) + np.min(self.datain)
		


		modelname = outdir+'model.fits'
		fits.writeto(outdir+'model.fits',model_image,self.head,overwrite=True)
		
		model_image[np.where(model_image==0)] = np.nan
#		if Norm == True:
#			model_image = model_image*(np.max(self.datain)-np.min(self.datain)) + np.min(self.datain)
		#set PA to constant

		

		residual = self.datain - model_image
		resname = outdir+'residuals.fits'
		fits.writeto(resname,residual,self.head,overwrite=True)

		return isolist


