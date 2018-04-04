__author__ = "Filippo Maccagni"
__copyright__ = "Fil8"
__email__ = "filippo.maccagni@gmail.com"

import sys, string, os
import numpy as np
import yaml
from kk import *
from radiobs import conv_units, cubeful, hi
from astropy import wcs
from astropy.io import fits, ascii
from astropy import units as u
from astropy.table import Table, Column, MaskedColumn
from astroquery.vizier import Vizier
from scipy import signal
import astropy.coordinates as coord
from mpdaf.obj import Spectrum, WaveCoord
from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib import rc


C=2.99792458e5 #km/s
HI=1.420405751e9 #Hz

####################################################################################################
#call class convert unit
cu = conv_units()
cubef = cubeful()
hi = hi()

class absex:
	'''
	
	Class for spectral studies (find continuum sources, extract spectra, analyze spectra)
	
	'''
	def __init__(self, file=None):
		'''
	
		Set logger for spectrum extraction
		Find config file
		If not specified by user load default.cfga
	
		'''
		
		cfg = open(file)
		self.cfg_par = yaml.load(cfg)

		import json

		key = 'general'

		self.workdir  = self.cfg_par[key].get('workdir', None)
		self.cubename = self.workdir + self.cfg_par[key].get('cubename', None)
		self.contname = self.workdir + self.cfg_par[key].get('contname', None)

		self.absdir  = self.workdir+'abs/'
		self.specdir = self.workdir+'spec/'
		self.plotdir = self.workdir+'plot/'

		print(json.dumps(self.cfg_par, indent=4, sort_keys=False))

		self.set_dirs()

	def enable_task(self,config,task):

		a = config.get(task, False)
		if a:
			return a['enable']
		else:
			False

	def set_dirs(self):
		'''
	 
		Sets directory strucure and filenames
		Creates directory abs/ in basedir+beam and subdirectories spec/ and plot/

		OUTPUT:
			tab : table of catalog
			flux : flux of continuum sources abouve set threshold
	 
		'''

		if os.path.exists(self.absdir) == False:
			 os.makedirs(self.absdir)                 
		if os.path.exists(self.specdir) == False:
			 os.makedirs(self.specdir)    
				
		if os.path.exists(self.plotdir) == False:
			 os.makedirs(self.plotdir)

	def source_catalog(self):
		key = 'source_catalog'
		width = self.cfg_par[key].get('width', '2')
		catalog = self.cfg_par[key].get('catalog', 'NVSS')
		thresh = self.cfg_par[key].get('thresh', 10e-3)
		centre = self.cfg_par[key].get('centre_coord',None)
		catalog_table = self.absdir + '/cat_src_absex.txt'

		Vizier.ROW_LIMIT = -1

		p = Vizier.query_region(coord.SkyCoord(centre[0],centre[1], unit=(u.hourangle, u.deg), 
			frame = 'icrs'), width = width, catalog = catalog)
		tab = p[0]

		ra_deg = []
		dec_deg = []
		if catalog == 'NVSS':
			for i in xrange (0, len(tab['RAJ2000'])):
			   tab['RAJ2000'][i] = string.join(string.split(tab['RAJ2000'][i],' '),':')
			   ra_deg.append(cu.ra2deg(tab['RAJ2000'][i]))
			   tab['DEJ2000'][i] = string.join(string.split(tab['DEJ2000'][i],' '),':')
			   dec_deg.append(cu.dec2deg(tab['DEJ2000'][i]))

			above_thresh = tab['S1.4']<thresh
		
		for i in xrange(1,len(tab.colnames)):
			tab[tab.colnames[i]][above_thresh] = np.nan

		tab =  Table(tab, masked=True)

		ascii.write(tab, catalog_table, overwrite=True)

		if self.cfg_par[key].get('simulate_continuum') == True:
			self.sim_cont_from_cube(catalog_table,catalog)

		return tab

	def abs_ex(self, verbose = True):

		cubefile = fits.open(self.cubename)  # read input
		hdr = cubefile[0].header
		sci = cubefile[0].data 
		sci = sci.squeeze()
		x = hdr['NAXIS1']
		y = hdr['NAXIS2']
		z = hdr['NAXIS3']
		cen_imx = hdr['CRPIX1']
		cen_imy = hdr['CRPIX2']
		freq0 = hdr['CRVAL3']
		freq_del = hdr['CDELT3']

		key = 'source_catalog'
		if self.enable_task(self.cfg_par,key):
			catalog_table = self.absdir + 'cat_src_absex.txt'
		   
			tab = ascii.read(catalog_table)
		
			if self.cfg_par[key].get('catalog', 'NVSS') == 'NVSS':

				self.J2000_name = tab['NVSS']
				self.ra = tab['RAJ2000']
				self.dec = tab['DEJ2000']
				self.flux_cont = tab['S1.4']*1e-3
				self.pixels = cu.coord_to_pix(self.cubename,self.ra,self.dec, verbose=False)
		
		key = 'abs_ex'

		src_id = np.arange(0,self.ra.size+1,1)
		freq = cubef.zaxis(self.cubename)

		self.abs_mean_rms = np.zeros(self.pixels.shape[0])
		self.abs_los_rms = np.zeros(self.pixels.shape[0])
		self.tau_los_rms = np.zeros(self.pixels.shape[0])
		outnames = []
		count_thresh =0
		count_fov = 0
		count_blanks = 0
		average_noise = []
		for i in xrange(0,self.pixels.shape[0]):

			# extract spectrum from each line of sight
			flux = np.zeros(freq.shape[0])
			madfm = np.zeros(freq.shape[0])

			if str(self.pixels[i,0]) == 'nan' or str(self.pixels[i,1]) == 'nan':
				count_thresh +=1
				pass

			elif (0 < int(self.pixels[i,0]) < x and
					0 < int(self.pixels[i,1]) < y): 
					
					pix_x_or = int(self.pixels[i,0])
					pix_y_or = int(self.pixels[i,1])
					for j in xrange(0, z):
						chrom_aber = self.cfg_par[key].get('chrom_aberration', False)
						#correct for chromatic aberration
						if chrom_aber == True:

							if (self.cfg_par[key].get('zunit','Hz') == 'm/s'):
								freq_real= freq* 1e2
								freq_real = (kk.C*kk.HI) /  (freq + kk.C)
								freq_real0 = (kk.C*kk.HI) /  (hdr['CRVAL3']*1e2 + kk.C)
								freq_del = (freq_real0 - freq_real[-1] )/ len(freq_real)
							#depending if the cube is in velocity or frequency ?
								scale = (freq_real0 - j*freq_del) / freq_real0
							pix_x = (pix_x_or - hdr['CRPIX1']) * scale + hdr['CRPIX1']
							pix_y = (pix_y_or - hdr['CRPIX1']) * scale + hdr['CRPIX1']
							pix_x = int(round(pix_x,0))
							pix_y = int(round(pix_y,0))
						else:
							pix_x = pix_x_or
							pix_y = pix_y_or
						
						if  (0 < pix_x < x and
							 0 < pix_y < y): 
							flux[j] = sci[j, pix_y, pix_x]
						else:
							flux[j] = 0.0
						
						# determine the noise of the spectrum [Whiting 2012 et al.] in each channel
						# MADMF: median absolute deviation from the median
						# extract a region were to determine the noise: A BOX around the l.o.s.
						if (pix_x+10 < hdr['NAXIS1'] and  pix_x-10 > 0 and
						   pix_y+10 < hdr['NAXIS2'] and pix_y - 10 > 0):
								rms = np.nanmedian(sci[j, pix_x +10:pix_x + 10, pix_y - 5:pix_y + 5])
								if rms != 0.0:
									med2 = np.abs(sci[j, pix_y, pix_x] - rms)
									madfm[j] = np.nanmedian(med2) / 0.6744888
								else:
									madfm[j] = 0.0
						else:
							madfm[j] = 0.0


						self.abs_mean_rms[i] = np.nanmean(madfm) 

					if np.nansum(flux) == 0.:
						count_blanks +=1
						if verbose == True:
							print '# Blank spectrum:\t'+str(src_id[i])+' '+self.J2000_name[i]+' #'
						continue

					# measure noise in the spectrum outside of the line
					end_spec = float(sci.shape[0])
					end_spec_th = int(end_spec/3.)
					end_spec = int(end_spec)
					mean_rms = (np.std(flux[0:end_spec_th]) +
										np.std(flux[end_spec-end_spec_th:end_spec])) / 2.
					mean_rms_arr = np.zeros(sci.shape[0])+mean_rms
					print self.flux_cont[i], mean_rms

					average_noise.append(mean_rms)

					tau = hi.optical_depth(flux,self.flux_cont[i])
					if np.nansum(madfm)!= 0.0:
						tau_noise = hi.optical_depth(madfm,self.flux_cont[i])
					else:
						tau_noise = np.zeros(sci.shape[0])

					#write spectrum
					out_spec = str(self.specdir+str(src_id[i])+'_J'+self.J2000_name[i])+'.txt'
					outnames.append(out_spec)

					flag_chans = self.cfg_par[key].get('flag_chans', None)
					if flag_chans != None:
						index_flags_l = (np.abs(freq - flag_chans[0])).argmin()
						for k in xrange(1,len(flag_chans)):
							index_flags = (np.abs(freq - flag_chans[k])).argmin()
							flux[index_flags_l:index_flags] = np.nan
							index_flags_l = index_flags

					if self.cfg_par[key].get('zunit','Hz') == 'm/s':
						xcol = 'Velocity [m/s]'
					elif self.cfg_par[key].get('zunit','Hz') == 'km/s':
						xcol = 'Velocity [km/s]'
					elif self.cfg_par[key].get('zunit','Hz') == 'MHz':
						xcol = 'Frequency [MHz]'
					else:
						xcol = 'Frequency [Hz]'

					t = Table([freq, flux, madfm, tau, tau_noise, mean_rms_arr], 
						names=(xcol,'Flux [Jy]','Noise [Jy]', 'Optical depth','Noise optical depth', 'Mean noise [Jy]'),
						meta={'name': 'Spectrum'})
					ascii.write(t,out_spec,overwrite=True)
					if verbose==True:
						print '# Extracted spectrum: \t' +str(src_id[i])+' '+self.J2000_name[i]+' #'

					polysub = self.cfg_par[key].get('polynomial_subtraction', False) 
					if polysub == True:

						deg = self.cfg_par[key].get('degree_pol',3)
						sub_flux = self.poly_sub(freq,flux,deg)
						sub_madfm = madfm.copy()
						sub_od = tau.copy()
						sub_noise_od = tau_noise.copy()

						out_spec_polysub = str(self.specdir+str(src_id[i])+'_J'+self.J2000_name[i])+'_psub.txt'
						t = Table([freq, sub_flux, sub_madfm, sub_od, sub_noise_od, mean_rms_arr], 
							names=(xcol,'Flux [Jy]','Noise [Jy]', 'Optical depth','Noise optical depth', 'Mean noise [Jy]'),
							meta={'name': 'Spectrum'})
						ascii.write(t,out_spec_polysub,overwrite=True)

					dohan = self.cfg_par[key].get('hanning', False)
					if dohan == True:
						
						window = self.cfg_par[key].get('window', 1)

						han_flux  = self.hanning_spec(flux)
						han_madfm = self.hanning_spec(madfm)
						han_od = self.hanning_spec(tau)
						han_noise_od = self.hanning_spec(tau_noise)
						
						out_spec_han = str(self.specdir+str(src_id[i])+'_J'+self.J2000_name[i])+'_han.txt'
						t = Table([freq, han_flux, han_madfm, han_od, han_noise_od, mean_rms_arr], 
							names=(xcol,'Flux [Jy]','Noise [Jy]', 'Optical depth','Noise optical depth', 'Mean noise [Jy]'),
							meta={'name': 'Spectrum'})
						ascii.write(t,out_spec_han,overwrite=True)

					if polysub == True and dohan == True:

						han_sub_flux  = self.hanning_spec(sub_flux)
						han_sub_madfm = self.hanning_spec(sub_madfm)
						han_sub_od = self.hanning_spec(sub_od)
						han_sub_noise_od = self.hanning_spec(sub_noise_od)
						
						out_spec_han = str(self.specdir+str(src_id[i])+'_J'+self.J2000_name[i])+'_psub_han.txt'
						t = Table([freq, han_sub_flux, han_sub_madfm, han_sub_od, han_sub_noise_od, mean_rms_arr], 
							names=(xcol,'Flux [Jy]','Noise [Jy]', 'Optical depth','Noise optical depth', 'Mean noise [Jy]'),
							meta={'name': 'Spectrum'})
						ascii.write(t,out_spec_han,overwrite=True)


		print '# Total number of sources: \t'+str(self.pixels.shape[0])
		print '# Sources flagged: \t\t'+str(count_thresh)
		print '# Blank spectra:\t\t'+str(count_blanks)
		print '# Total number of spectra: \t'+str(self.pixels.shape[0]-count_thresh-count_fov-count_blanks)
		print '# Average noise in spectra: \t'+str(round(np.nanmean(average_noise)*1e3,1))+' mJy/beam'

		return 0 

	def abs_plot(self, spec_src_name):
		'''		
		Plots spectra of all radio sources found by find_src_imsad 
		saved in basedir/beam/abs/spec.
		Plots are stored in basedir/beam/abs/plot
		
		IN
			Spectra extracted by spec_ex
		
		IN cfga
			abs_ex_plot_xaxis= ' '      #: X-axis units ['velocity','frequency'] 
			abs_ex_plot_yaxis= ' '      #: Y axis units ['flux','optical depth']
			abs_ex_plot_redsrc= True    #: plots line at redshift of source in spectrum redshift must be stored in table of load_src_csv
			abs_ex_plot_title= True     #: plot title: J2000 name of radio source
			abs_ex_plot_format= ' '     #: format of plot ['.pdf','.jpeg','.png']
		
		OUT
			For each source outputs have the following name:
			J2000_xaxis-unit_yaxis-unit.plot_format = J220919.87+180920.17_vel_flux.pdf

		'''

		key = 'abs_plot'

		os.chdir(self.specdir)
		
		params = {
				  'text.usetex': True,
				  'text.latex.unicode': True
				   }
		rc('font', **{'family': 'serif', 'serif': ['serif']})        
		plt.rcParams.update(params)
		
		for i in xrange(0,len(np.atleast_1d(spec_src_name))):
			
			#load data and labels 
			spec_name = spec_src_name[i]
			if os.path.isfile(spec_name) == True:
			
				# Set plot specs
				font_size = 16            
				plt.ioff()
				fig = plt.figure(figsize =(8,6))
				fig.subplots_adjust(hspace=0.0)
				gs = gridspec.GridSpec(1, 1)
				plt.rc('xtick', labelsize=font_size-2)
				plt.rc('ytick', labelsize=font_size-2) 

				# Initialize subplots
				ax1 = fig.add_subplot(gs[0])
				ax1.set_xlabel('')
				ax1.set_ylabel('') 
				
			
				spec_vec = ascii.read(spec_name)


				x_data = np.array(spec_vec[spec_vec.colnames[0]],dtype=float)

				flag_chans = self.cfg_par['abs_ex'].get('flag_chans', None)
				if self.cfg_par['abs_ex'].get('zunit') == 'm/s':
					x_data /= 1e3

				y_data = np.array(spec_vec[spec_vec.colnames[1]],dtype=float)*1e3
				y_sigma = np.array(spec_vec[spec_vec.colnames[2]])

				ax1.set_xlabel(r'$cz\,(\mathrm{km}\,\mathrm{s}^{-1})$', fontsize=font_size)
				ylabh = ax1.set_ylabel(r'S\,$[\mathrm{mJy}\,\mathrm{beam}^{-1}]$', fontsize=font_size)
				ylabh.set_verticalalignment('center')

				# Calculate axis limits and aspect ratio
				x_min = np.min(x_data)
				x_max = np.max(x_data)
				y1_array = y_data[np.where((x_data>x_min) & (x_data<x_max))]
				y1_min = np.nanmin(y_data)*1.1
				y1_max = np.nanmax(y_data)*1.1

				# Set axis limits
				ax1.set_xlim(x_min, x_max)
				ax1.set_ylim(y1_min, y1_max)
				ax1.xaxis.labelpad = 6
				ax1.yaxis.labelpad = 10

				# Plot spectra 
#                if self.abs_ex_plot_linestyle == 'step':
					#ax1.plot(x_data, y_data, color='black', linestyle='-')

				ax1.step(x_data, y_data, where='mid', color='black', linestyle='-')
				


				# Plot noise
				ax1.fill_between(x_data, -y_sigma, y_sigma, facecolor='grey', alpha=0.5)

				if flag_chans != None:
					flag_chans = np.array(flag_chans) 
					if 	self.cfg_par['abs_ex'].get('zunit') == 'm/s':
						flag_chans = np.divide(flag_chans,1e3)
					index_flags_l = (np.abs(x_data - flag_chans[0])).argmin()
					for k in xrange(1,len(flag_chans)):
						index_flags = (np.abs(x_data - flag_chans[k])).argmin()
						ax1.fill_between([x_data[index_flags_l],x_data[index_flags]],y1_min,y1_max, facecolor='grey', alpha=0.3)



				# Plot stuff
				ax1.axhline(color='k', linestyle=':', zorder=0)
				
				redshifts = self.cfg_par[key].get('redshift_sources',None)
				if len(redshifts) == 2:
						ax1.fill_between([redshifts[0],redshifts[1]],y1_min,y1_max, facecolor='red', alpha=0.1)

				# Add title        
				#if self.abs_ex_plot_title == True:
				#	ax1.set_title('%s' % (self.J2000_name[i]), fontsize=font_size+2) 
				#ax1.axes.titlepad = 8

				# Add minor tick marks
				ax1.minorticks_on()

				# Save figure to file
				outplot = os.path.basename(spec_src_name[i])
				outplot = string.split(outplot,'.')[0]
				outplot = self.plotdir+outplot+'.png'
				plt.savefig(outplot,
							overwrite = True)
				print '# Plotted spectrum of source ' + os.path.basename(spec_src_name[i])+'. #'
			else:
				print '# Missing spectrum of source ' + os.path.basename(spec_src_name[i])+'. #'

# #	def hanning_spec(self,freq,flux, window):
# 		"""smooth the data using a window with hanning window of requested size.
		
# 		This method is based on the convolution of a scaled window with the signal.

# 		input:
# 			x: the input signal 
# 			window_len: the dimension of the smoothing window; should be an odd integer
# 			window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
# 				flat window will produce a moving average smoothing.

# 		output:
# 			the smoothed signal

# 		TODO: the window parameter could be the window itself if an array instead of a string
# 		NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
# 		"""

# 		# if flux.ndim != 1:
# 		# 	raise ValueError, "smooth only accepts 1 dimension arrays."

# 		# if flux.size < window:
# 		# 	raise ValueError, "Input vector needs to be bigger than window size."


# 		# if window<3:
# 		# 	return flux

# 		# old_len = len(flux)

# 		# s=np.r_[flux[window-1:0:-1],flux,flux[-2:-window-1:-1]]
		
# 		# w =np.hanning(window)

# 		# flux=np.convolve(w/np.sum(w),s, mode='valid')

# 		# flux = flux[(window/2-1):-(window/2)]	

# 		# new_len = len(flux)


# 		# if (old_len < new_len): 
# 		# 	diff = new_len - old_len
# 		# 	flux =  flux[0:-diff]
# 		# elif old_len > new_len:
# 		# 	diff = old_len - new_len
# 		# 	listdiff = np.empty([diff])
# 		# 	listdiff[:] = np.nan
# 		# 	flux = np.append(flux,listdiff)

# 		# step_freq = (freq[-1]-freq[0])/len(freq)
# 		# new_step_freq = step_freq*window
# 		# new_freq = np.arange(freq[0],freq[-1]+new_step_freq,new_step_freq)
# 		# flux = np.interp(new_freq, freq, flux)



# 		return  new_freq, flux

	def hanning_spec(self,flux):

		new_flux = flux.copy()
		new_flux[0] = (flux[0]+flux[1])/2.
		for i in xrange(1,len(flux)-1):

			new_flux[i] = (flux[i-1]+2.*flux[i]+flux[i+1])/4.

		new_flux[-1] = (flux[-2]+flux[-1])/2.

		return new_flux

	def poly_sub(self,x,y,degree):

		if self.cfg_par['abs_ex'].get('zunit') == 'm/s':
			unit_z = u.m / u.s

		step = (x[-1]-x[0])/len(x)
		wave_for_spec = WaveCoord(cdelt=step, crval=x[0], cunit= unit_z)
		spe = Spectrum(wave=wave_for_spec, data=y)
		cont = spe.poly_spec(degree)

		cont_sub = y - cont.data

		return cont_sub

	def sim_cont_from_cube(self, table,catalog):

		tab = ascii.read(table)

		if catalog == 'NVSS':

				ra = tab['RAJ2000']
				dec = tab['DEJ2000']
		
				major = tab['MajAxis']
				minor = tab['MinAxis']		

				angle1=np.radians(0.0)
				cosangle1=np.cos(angle1)
				sinangle1=np.sin(angle1)

		pixels = cu.coord_to_pix(self.cubename,ra,dec,verbose=False)


		cubefile = fits.open(self.cubename)
		contdata = cubefile[0].data
		cubehead = cubefile[0].header

		if cubehead['NAXIS'] > 3:
			contdata = np.zeros([contdata.shape[2],contdata.shape[3]])

			del cubehead['CTYPE4']
			del cubehead['CDELT4']    
			del cubehead['CRVAL4']
			del cubehead['CRPIX4']
			del cubehead['NAXIS4']
		else:
			contdata = np.zeros([contdata.shape[1],contdata.shape[2]])

		del cubehead['CRPIX3'] 
		del cubehead['CRVAL3']
		del cubehead['CDELT3']
		del cubehead['CTYPE3']
		del cubehead['NAXIS3']
		del cubehead['NAXIS']

		w=wcs.WCS(cubehead)    

		xnum = np.linspace(0,cubehead['NAXIS1'],cubehead['NAXIS1'])
		ynum = np.linspace(0,cubehead['NAXIS2'],cubehead['NAXIS2'])	
		x, y = 	np.meshgrid(xnum, ynum)

		for i in xrange(0,pixels.shape[0]):

			xc=float(pixels[i][0])
			yc=float(pixels[i][1])

			if str(xc) == 'nan' or str(yc) == 'nan':
				pass
			else:
				if minor[i]/3600. >= float(cubehead['CDELT2']) and major[i]/3600. >= float(cubehead['CDELT2']):
					a = major[i]/3600./float(cubehead['CDELT2'])/2.
					b = minor[i]/3600./float(cubehead['CDELT2'])/2.

					ell = np.power(x-xc, 2)/np.power(a,2) + np.power(y-yc, 2)/np.power(b,2)
					index_ell = np.where(np.less_equal(ell,1))
					contdata[index_ell] = 1
				else:
					contdata[int(yc),int(xc)] = 1

		fits.writeto(self.absdir+'tmp_cont.fits',contdata,cubehead,overwrite=True)






