__author__ = "Filippo Maccagni"
__copyright__ = "Fil8"
__email__ = "filippo.maccagni@gmail.com"

import sys,string,os,math
import numpy as np
import Cosmo as c
from kk import *
import aipy
import ephem, datetime
from astropy import units as u
from astropy.coordinates import Angle, AltAz, SkyCoord, EarthLocation
from astropy.time import Time, TimeDelta

class observe:
	'''
	+----+
	  observe
	+----+
		Tools providing info on radio observations
	 	
	 	- pbeam_FWHM
	 		full-width at half maximum iof the primary beam of a MeerKAT observation
	 	- theo_rms
	 		theoretical r.m.s. for a MeerKAT or Apertif observation, 
	 		given the number of antennas and integration time
	 	- local_standard_time
	 		LST of a Westerbork observation given the UCT time of the observation
	 	- alt_az 
	 		ALT/AZ of an observed source in a MIRIAD format observation
	'''


	def __init__(self):

		#define constants used in this module
		self.C= kk.C

	def pbeam_FWHM(self,telescope, obs_freq):

		if telescope == 'meerkat' or telescope == 'MeerKAT' or telescope == 'meerKAT' or telescope == 'meer':
			# dish diameter (m)
			diam = 13.5

		# Primary beam FWHM (deg)
		pbeam_fwhm = 1.02*(self.C*1e-2)/obs_freq/diam/np.pi*180.

		return pbeam_fwhm

	def theo_rms(self,telescope,nant,bw,int_time):
		'''
		Determines the predicted noise per frequency channel for the perfomance of a given telescope
		int_time is in seconds
		bw is in ?????
		'''
		
		if telescope == 'meerkat' or telescope == 'MeerKAT' or telescope == 'meerKAT' or telescope == 'meer':
			coreta = 0.88
			tsys = 100.
			jyperk = 2.
		elif telescope == 'apertif' or telescope == 'Apertif' or telescope == 'APERTIF' or telescope == 'ape':
			coreta = 0.88
			tsys = 100.
			jyperk = 2.
		else:
			print 'Telescope not known'
			return 0

		tsys = tsys*coreta
		noise_freq = tsys*jyperk/np.sqrt(nant*(nant-1)*np.abs(bw)*int_time)

		noise_freq *= 1e3 #mJy/beam

		print 'Expected noise per channel = '+str(round(noise_freq,3))+' mJy/beam'

		return noise_freq

	def local_standard_time(self, telescope, timeobs):

		sun=ephem.Sun()
		# Define the location of APERTIF
		if telescope == 'apertif' or telescope == 'Apertif' or telescope == 'APERTIF' or telescope == 'ape' or telescope == 'WSRT' or telescope == 'Westerbork':
			loc = ephem.Observer()
			WSRT_ra_deg= Angle('6.60334d')
			WSRT_ra = WSRT_ra_deg.to_string(unit=u.degree, sep=':')
			loc.lon = WSRT_ra
			WSRT_dec_deg= Angle('52.91474d')
			WSRT_dec = WSRT_dec_deg.to_string(unit=u.degree, sep=':')
			loc.lon = WSRT_dec
			loc.elevation = 17
			loc.horizon = '10:00:00' # degrees

		loc.date = timeobs.datetime
		date = str(timeobs.datetime.year)+'-'+str(timeobs.datetime.month)+'-'+str(timeobs.datetime.day)
		sun.compute(loc)
		# sidereal time == ra (right ascension) is the highest point (noon)
		lst = loc.sidereal_time() - sun.ra
		lst = ephem.hours(lst + ephem.hours('12:00')).norm
		lst = str(lst)
		lst = string.split(lst,':')
		lst_full = date+'T'+str(int(lst[0]))+':'+str(int(lst[1]))+':'+str(lst[2])
		j2000_jd = 2451544.5
		UT_dec = float(timeobs.datetime.hour +timeobs.datetime.minute/60.+ timeobs.datetime.second/3600.)
		LST = 100.46 + 0.985647 * (float(timeobs.jd)-j2000_jd) + WSRT_ra_deg.degree + 15.*(UT_dec) / 360.
		#error of 15 s within 100 years from J2000
		LST = Angle(LST, unit=u.deg)

		print LST.to_string(unit=u.hour, sep=':')

		return lst_full

	def alt_az(self,telescope,observation):

		#open miriad observation
		uv = aipy.miriad.UV(str(observation))

		#count visibilities
		ant1=0
		ant2=1
		uvnum=0
		for preamble, data, flags in uv.all(raw=True):
			#    #count number of visibilities per baseline
					if ( (preamble[2][0] == ant1) and (preamble[2][1]==ant2)) :
						uvnum+=1


		#read Ra & Dec of observation
		ra = Angle(uv['obsra'], u.radian)
		dec = Angle(uv['obsdec'], u.radian)
		coord = SkyCoord(ra.degree, dec.degree,  unit=(u.deg, u.deg))
		#time and date of observation
		end_time=Time(uv['time'], format = 'jd', scale='utc')
		print end_time.datetime
		inttime=uv['inttime']
		inttime = TimeDelta(inttime*uvnum/4.01, format='sec')
		start_time = end_time-inttime
		print start_time.jd, end_time.jd
		#define location
		if telescope == 'apertif' or telescope == 'Apertif' or telescope == 'APERTIF' or telescope == 'ape' or telescope == 'WSRT' or telescope == 'Westerbork':

			WSRT_ra_deg= Angle('6.60334d')
			WSRT_dec_deg= Angle('52.91474d')
			#define alternative location
			WSRT = EarthLocation(lat=WSRT_ra_deg, lon=WSRT_ra_deg, height=17*u.m)

		else:
			print 'Telescope not known'
			return 0


		#define alternative location
		WSRT_ra_deg= Angle('6.60334d')
		WSRT_dec_deg= Angle('52.91474d')
		WSRT = EarthLocation(lat=WSRT_ra_deg, lon=WSRT_ra_deg, height=17*u.m)
		j2000_jd = 2451544.5

		UT_decS = float(start_time.datetime.hour +start_time.datetime.minute/60.+ start_time.datetime.second/3600.)
		LSTS = (100.46 + 0.985647 * (float(start_time.jd)-j2000_jd) + WSRT_ra_deg.degree + 15.*(UT_decS) )/ 360.
		#error of 15 s within 100 years from J2000
		LSTS = Angle(LSTS, unit=u.deg)
		LSTS = LSTS.to_string(unit=u.hour, sep=':')

		UT_decE = float(end_time.datetime.hour +end_time.datetime.minute/60.+ end_time.datetime.second/3600.)
		print UT_decE, end_time.datetime.hour
		LSTE = (100.46 + 0.985647 * (float(end_time.jd)-j2000_jd) + WSRT_ra_deg.degree + 15.*(UT_decE) )/ 360.
		#error of 15 s within 100 years from J2000
		LSTE = Angle(LSTE, unit=u.deg)
		LSTE = LSTE.to_string(unit=u.hour, sep=':')
		print UT_decS, UT_decE
		print LSTS, LSTE
		#Local Standard Time
		#s_lst = self.local_standard_time('WSRT',start_time)
		#e_lst = self.local_standard_time('WSRT',end_time)
		times = [LSTS,LSTE]
		#starting frame
		frame = AltAz(obstime=times, location=WSRT)

		obs_altasz = coord.transform_to(frame)

		#hour angle 
		#e_lha = string.split(e_lst,'T')
		#e_lha = Angle(e_lha[1], u.hourangle)
		#e_lha = e_lha-ra

		return obs_altasz
