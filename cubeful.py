__author__ = "Filippo Maccagni"
__copyright__ = "Fil8"
__email__ = "filippo.maccagni@gmail.com"

import sys, string, os
import numpy as np
import pyregion
from kk import *
from astropy.io import fits, ascii
from astropy import units as u
from astropy import wcs
from astropy.coordinates import Angle
from astropy import units as u
from astropy.table import Table

class cubeful:
	'''
	+----+
	  cubeful
	+----+
	  Tools to handle cube data and headers
		- zcut 
			divide cube in subchuncks along the z axis
	'''

	def __init__(self):

		self.C = kk.C
		self.HI = kk.HI

	def zcut(self, cube, width):
		'''
		Divides cube in subchuncks along the z-axis of width specified by user
		INPUT
			cube: cube to divide in .fits format
			width: width of each subchunkcs (must be of same unit as z-axis) 
		OUTPUT
			cube_startfreq-endfreq.fits: subchuncks of initial cube
		'''

		#save masked continuum image 
		indir = os.path.dirname(cube)
		outdir = indir+'/cubelets/'
		if os.path.exists(outdir) == False:
			os.mkdir(outdir)

		cubein = fits.open(cube)
		dc = cubein[0].data
		dh = cubein[0].header

		freq = (np.linspace(1, dc.shape[0], dc.shape[0]) - dh['CRPIX3']) * dh['CDELT3'] + dh['CRVAL3']

		start_dc = freq[0]
		end_dc = freq[-1]
		zunit = dh['CTYPE3']
		del dh['FRAME']
		width = string.split(width,' ')
		wunit = str(width[1])
		width = float(width[0])
		width *= 1e6

		zint = np.arange(start_dc,end_dc+width,width)
		indint =[]
		for i in xrange (0,len(zint)):
			indint.append(np.abs(freq - zint[i]).argmin())

		zint /= 1e6
		unit = 'MHz'
		
		for i in xrange(1,len(indint)):
			sub_dc = dc[indint[i-1]:indint[i],:,:]
			sub_name = outdir+'cube_'+str(int(zint[i-1]))+'-'+str(int(zint[i]))+unit+'.fits'

			sub_dh = dh.copy()
			sub_dh['CRVAL3'] = zint[i-1]
		
			sub_dh['CRVAL3'] *= 1e6
		
			sub_dh['NAXIS3'] = indint[i]-indint[i-1]
			fits.writeto(sub_name,sub_dc,sub_dh,overwrite=True)

		print '''\t+---------+\n\t Cubes cut\n\t+---------+'''


	def zsum(self,cube):
		'''
		Sums cube over z axis to create a continuum image / mom 0 map
		INPUT
			cube 
		OUTPUT
			cont_im_midfreq.fits: continuum image of the initial cube
		'''		

		#save masked continuum image 
		indir = os.path.dirname(cube)
		outdir = indir+'/../cont_ims/'
		if os.path.exists(outdir) == False:
			os.mkdir(outdir)

		cubein = fits.open(cube)
		dc = cubein[0].data
		dh = cubein[0].header

		zunit = dh['CTYPE3']
		if zunit == 'Hz':
			delta = float(dh['CDELT3'])/1e6
			start = float(dh['CRVAL3'])/1e6

		width = float(dc.shape[0])
		mid_freq = start+delta*width/2.

		cont_im = outdir+'cont_'+str(round(mid_freq,0))+'.fits'
		dh['NAXIS']= 2
		del dh['NAXIS3']
		del dh['CTYPE3']
		del dh['CRVAL3']
		del dh['CRPIX3']
		del dh['CDELT3']
		del dh['CUNIT3']

		index = np.isnan(dc)
		dc[index] = 0.0

		ic = np.zeros([dc.shape[1],dc.shape[2]])
		nan_counter = 0
		for i in xrange(0,dc.shape[0]):

			ic[:,:] += dc[i,:,:]
			if np.nanmean(dc[i,:,:]) == 0:
				nan_counter +=1
		#for i in xrange(0,dc.shape[2]):
		#	for j in xrange(0,dc.shape[1]):
		#		ic[j,i] = np.nanmean(dc[:,j,i])
		
		#ic = np.nanmean(dc,axis=0)
		#print nan_counter, os.path.basename(cube)
		ic = np.divide(ic,width-nan_counter)
		print nan_counter
		fits.writeto(cont_im,ic,dh,overwrite=True)

	def zaxis(self,cubename):

		cubefile = fits.open(cubename)  # read input
		hdr = cubefile[0].header
				
		freq = (np.linspace(1, hdr['NAXIS3'], hdr['NAXIS3']) - hdr['CRPIX3']) * hdr['CDELT3'] + hdr['CRVAL3']
 
		return freq           

	#def merge_to_cube()

	def zcompress(self,cubelist,outcube):
		'''
		Merge continuum images in one single cube
		INPUT
			cubelist : list of paths of continuum images to merge
		OUTPUT
			in final directory of cubelist
			outcube: full path output cube]
			outtable: table with major and minor axis of PSF/dirty beam of continuum images
					  (only if PSF are different from eachother)
		
		RETURN
			table: table of beams
			
		'''				

		#size of cube, since all continuum images may not have same NAXIS1/2 lenght
		#but must have same pixel size, otherwise use racont.reproj_image
		lenaxis=[]
		for i in xrange(0,len(cubelist)):
			base = fits.open(cubelist[0])
			basedata = base[0].data
			baseheader =  base[0].header
			lenaxis.append(float(baseheader['NAXIS2']))
		

		cubeshape = np.min(lenaxis)
		print cubeshape

		cube = fits.PrimaryHDU()

		cube.data = np.zeros([len(cubelist),cubeshape,cubeshape])

		bmin = []
		bmaj = []
		cubenames = []
		print cube.data.shape
		for i in xrange(0,len(cubelist)):
			base_tmp = fits.open(cubelist[i])[0]
			print base_tmp.data.shape

			y_diff = base_tmp.data.shape[0] - cube.data.shape[1]
			x_diff = base_tmp.data.shape[1] - cube.data.shape[2]
	
			y_appr =  y_diff/2 + y_diff/2		
			x_appr =  x_diff/2 + x_diff/2		
			print y_diff,y_appr,x_diff,x_appr
			if y_diff == y_appr and x_diff == x_appr:	
				cube.data[i,:,:] = base_tmp.data[y_diff/2:base_tmp.data.shape[0]-y_diff/2,
											x_diff/2:base_tmp.data.shape[1]-x_diff/2]
			else:
				y_newdiff = y_diff-y_appr
				x_newdiff = x_diff-x_appr
				cube.data[i,:,:] = base_tmp.data[y_diff/2:base_tmp.data.shape[0]-y_diff/2-y_newdiff,
											x_diff/2:base_tmp.data.shape[1]-x_diff/2-x_newdiff]



			cubename = string.split(cubelist[i],'/')
			cubenames.append(cubename[-1])
			bmin.append(base_tmp.header['BMIN'])
			bmaj.append(base_tmp.header['BMAJ'])

		cube.header['CRPIX3'] = 1
		cube.header['CDELT3'] = 1
		cube.header['CRVAL3'] = 1
		cube.header['CRPIX1'] = baseheader['CRPIX2']
		cube.header['CDELT1'] = baseheader['CDELT2']
		cube.header['CRVAL1'] = baseheader['CRVAL2']		
		cube.header['CTYPE1'] = baseheader['CTYPE1']		
		cube.header['CRPIX2'] = baseheader['CRPIX2']
		cube.header['CDELT2'] = baseheader['CDELT2']
		cube.header['CRVAL2'] = baseheader['CRVAL2']
		cube.header['CTYPE2'] = baseheader['CTYPE2']		
		cube.header['CTYPE3'] = 'CHAN'

		fits.writeto(outcube,cube.data,cube.header,overwrite=True)
		
		outtablenames = string.split(outcube,'.fits')
		outtable = outtablenames[0]+'_beams.txt'
		ascii.write([cubenames, bmin, bmaj], outtable, names=['ID', 'bmin_deg', 'bmaj_deg'])
		t = Table([cubenames, bmin, bmaj], names = ['ID', 'bmin_deg', 'bmaj_deg'] )


		return t












