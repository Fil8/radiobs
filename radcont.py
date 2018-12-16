__author__ = "Filippo Maccagni"
__copyright__ = "Fil8"
__email__ = "filippo.maccagni@gmail.com"

import sys, string, os
import numpy as np
import pyregion
from kk import *
from astropy.io import fits
from astropy import units as u
from astropy import wcs
from reproject import reproject_interp, reproject_exact, reproject_from_healpix
from astropy.utils.data import get_pkg_data_filename
from astropy.coordinates import Angle
from astropy import units as u
import montage_wrapper as montage

class radcont:
	'''
	+----+
	  radcont
	+----+
	  Tools for radio continuum observations
		- int_flux 
			measure integral flux of a source from a continuum image in .fits format
			within a region drawn by hand (using ds9)
		- mask_cont_forflux
			create masked continuum image: 0-values outside cutoff set by user
	'''

	def __init__(self):

		self.C = kk.C
		self.HI = kk.HI

	def flux_in_mask(self,maskdata,maskhead,cutoff,region):

		# set polygonal mask from ds9 region
		r = pyregion.open(region).as_imagecoord(maskhead)
		shape = (maskhead['NAXIS2'], maskhead['NAXIS1'])
		m = r.get_mask(shape=shape)
		
		background = np.nanmean(maskdata[m==False])
		noise = np.nanstd(maskdata[m==False])


		maskdata[m==False] = np.nan


		maskdata[maskdata<cutoff] = np.nan

		mask = maskdata.copy()
		mask[~np.isnan(mask)] = 1

		noise = np.multiply(noise,np.sqrt(np.sum(mask)))
		
		#maskdata -= background


		return maskdata, noise, background, np.nansum(mask)


	def int_flux(self, cont_im, cutoff, region, pbeam_corr=True):
		'''
		Measures the total flux density of a source in a continuum image 
		within a region drawn by hand (using ds9)
		INPUT
			cont_im: continuum image in .fits format
			cutoff: cutoff level above which measure the flux [units of image (Jy/beam)] 
		OUTPUT
			int_flux: total flux density of the selected source [Jy]
			cont_im_flux.fits: image of cutout where flux was measured (stored in workdir+int_flux_output/)
		'''


		cont_imf = fits.open(cont_im)
		conthead = cont_imf[0].header
		contdata = cont_imf[0].data

		maskdata = contdata.copy()
		mask_cont_data = self.mask_cont_forflux(maskdata, conthead, cutoff, region)

		contdata[mask_cont_data!=1.0] = np.nan

		#save masked continuum image 
		inputdir = os.path.dirname(cont_im)
		outputdir = inputdir+'/../int_flux_output/'
		if os.path.exists(outputdir) == False:
			os.mkdir(outputdir)

		outname = os.path.basename(cont_im)
		outname = string.strip(outname)
		outname = string.split(outname,'.fits')
		outfile = outputdir+outname[0]+'_flux.fits'
		fits.writeto(outfile,contdata,conthead,overwrite=True)

		#measure flux density
		total_flux = np.nansum(contdata)
		
		#correct for primary beam
		pix_area = -float(conthead['CDELT1'])*float(conthead['CDELT2'])

		if pbeam_corr==True:
			bmaj = float(conthead['BMAJ'])
			bmin = float(conthead['BMIN'])
			
			bmaj = Angle(bmaj,unit=u.degree)
			bmin = Angle(bmin,unit=u.degree)
			
			beam_area = 2.*np.pi*bmaj.degree/2.35482*bmin.degree/2.35482
			number_pix_beam= beam_area/pix_area

			total_flux /= number_pix_beam 

		total_flux = total_flux

		print str(os.path.basename(cont_im))+':\t total flux density = '+str(total_flux)

		return total_flux, pix_area



	def mask_cont_forflux(self, maskdata, maskhead, cutoff, region):

		# set polygonal mask from ds9 region
		r = pyregion.open(region).as_imagecoord(maskhead)
		shape = (maskhead['NAXIS1'], maskhead['NAXIS2'])
		print shape
		m = r.get_mask(shape=shape)
		maskdata[m==False] = np.nan

		maskdata[maskdata<cutoff] = 0.0
		maskdata[maskdata>=cutoff] = 1.0		
		
		return maskdata


	def read_psf(self, filename):
		'''
		Read psf dimensions from BMAJ, BMIN keywords of fits file
		INPUT
			filename : name of .fits image/cube
		OUTPUT
			[filename,bmin,bmaj,cd1,cd2,naxis1,naxis2]: numpy array with filename, major/minor axis 
														of psf, pixel_size, shape of .fits
			bmaj and bmin are given in arcminutes
		'''		
		filelist = fits.open(filename)
		filehead = filelist[0].header
		filedata = filelist[0].data
		
		bmaj = 0
		bmin = 0
		cd1 = 0
		cd2 = 0
		naxis1 = 0
		naxis2 = 0

		if 'BMAJ' in filehead:
			bmaj = float(filehead['BMAJ'])
		if 'BMIN' in filehead:
			bmin = float(filehead['BMIN'])
		if 'CDELT1' in filehead:
			cd1 = float(filehead['CDELT1'])
		elif 'CD1_1' in filehead:
			cd1 = float(filehead['CD1_1'])
		if 'CDELT2' in filehead:
			cd2 = float(filehead['CDELT2'])
		elif 'CD2_2' in filehead:
			cd2 = float(filehead['CD2_2'])

		naxis1 = float(filehead['NAXIS1'])
		naxis2 = float(filehead['NAXIS2'])

		cd1 = Angle(cd1, unit=u.degree)
		cd2 = Angle(cd2, unit=u.degree)
		bmaj = Angle(bmaj,unit=u.degree)
		bmin = Angle(bmin,unit=u.degree)

		array_tmp = [bmaj.arcminute, bmin.arcminute, cd1.arcsecond, cd2.arcsecond, naxis1, naxis2]

		return array_tmp


	def reproj_image(self,basename,slavename,outname,fluxtype):
		'''
		Merge continuum images in one single cube
		INPUT
			basename : list of paths of continuum images to merge
			slavename : list of paths of continuum images to merge			
			outname: 
			fluxtype: 
		OUTPUT
			in final directory of cubelist
			outcube:  full path output cube]
			outtable: table with major and minor axis of PSF/dirty beam of continuum image
					  (only if PSF are different from eachother)
		
		RETURN
			cubenames: list filenames of images that have been concatenated
			
		'''	


		base = fits.open(basename)[0]
		slave = fits.open(slavename)[0]
		if float(base.header['NAXIS']) == 4:
			base.header['NAXIS'] = 2
			del base.header['NAXIS3']
			del base.header['NAXIS4']
			del base.header['CRVAL3']
			del base.header['CDELT3']
			del base.header['CRPIX3']
			del base.header['CTYPE3']
			del base.header['CRVAL4']
			del base.header['CDELT4']
			del base.header['CRPIX4']
			del base.header['CTYPE4']
			base.data = np.squeeze(base.data)
			base.data = np.squeeze(base.data)
		elif float(base.header['NAXIS']) == 3:
			base.header['NAXIS'] = 2
			del base.header['NAXIS3']
			del base.header['CRVAL3']
			del base.header['CDELT3']
			del base.header['CRPIX3']
			del base.header['CTYPE3']
			base.data = np.squeeze(base.data)
		
		if 'WCSAXES' in slave.header:
			slave.header['WCSAXES'] = 2
			del slave.header['PC3_1']
			del slave.header['PC3_2']
			del slave.header['PC1_3']
			del slave.header['PC2_3']
			del slave.header['PC3_3']
		if 'CROTA3' in base.header:
			del base.header['CROTA3']
		if 'CROTA4' in base.header:
			del base.header['CROTA4']
		
		# if 'CD1_1' in slave.header:
		# 	pix_area = -float(slave.header['CD1_1'])*float(slave.header['CD2_2'])
		
		# elif 'CDELT1' in slave.header:
		# 	pix_area = -float(slave.header['CDELT1'])*float(slave.header['CDELT2'])

		# if 'BMAJ' in slave.header:
		# 	bmaj = float(slave.header['BMAJ'])
		# 	bmin = float(slave.header['BMIN'])
		# else:
		# 	bmaj = 0.0
		# 	bmin = 0.0
		# bmaj = Angle(bmaj,unit=u.degree)
		# bmin = Angle(bmin,unit=u.degree)
		

		if fluxtype == 'exact':
			newslave, footprint = reproject_exact(slave, base.header)
			#base.header['BMAJ'] = slave.header['BMAJ']
			#base.header['BMIN'] = slave.header['BMIN']
			fits.writeto(outname, newslave, base.header, clobber=True)
		if fluxtype == 'healpix':
			print slavename
			#filename_in = get_pkg_data_filename(slavename)
			filename_in = slavename
			newslave, footprint = reproject_from_healpix(filename_in, base.header)
			fits.writeto(outname, newslave, base.header, clobber=True)
		elif fluxtype == 'fast':
			newslave, footprint = reproject_interp(slave, base.header)
			base.header['BMAJ'] = slave.header['BMAJ']
			base.header['BMIN'] = slave.header['BMIN']
			fits.writeto(outname, newslave, base.header, clobber=True)
		elif fluxtype == 'montage':
			headername = 'tmp.hdr'
			montage.mGetHdr(basename,headername)
			montage.mProject(slavename,outname,headername)
			os.remove(headername)
		elif fluxtype == 'montage_healpix':
			fitsname = 'tmp.fits'
			headername = 'tmp.hdr'
			a = montage.mGetHdr(basename,headername)
			command = 'HPXcvt '+slavename+' '+fitsname
			os.system(command)
			#montage.HPXcvt(slavename,fitsname)
			command = 'mProject '+fitsname+' '+outname+' '+headername
			#montage.mProject(slavename,outname,headername)
			os.system(command)
			os.remove(fitsname)
			os.remove(headername)


		return 0







