__author__ = "Filippo Maccagni"
__copyright__ = "Fil8"
__email__ = "filippo.maccagni@gmail.com"

import sys,string,os,math
import numpy as np
import Cosmo as c
from astropy import units
from astropy.io import fits

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

KnhiA = 1.8216E18
KnhiE = 3.1E17
M26 = -23.33
STULLY = -9.64



class hi:

	def __init__(self):
		self.hi=HI
		self.m26=M26
		self.s_tully= STULLY
		self.knhi= Knhi
		self.nhiem= KnhiE
		self.T=TSPIN
		self.mhi=MHI
		self.msun=MSUN
		self.chi=CHI
		self.mp=MP

	def hiline(self,z):
		'''
		Estimates the expected velocity of the HI line given the redshift
		of the source
		INPUT
			z: redshift (float)
		OUTPUT
			hi.velhi: velocity in km/s
			hi.freq: frequency of the HI line in MHz
		'''

		freq=self.hi/(1+z)/1e06	#MHz
		velocity=C*((self.hi-self.freq)/self.freq)/1e5	#km/s

		print 'HI expected frequency = '+str(round(freq,3))+' mJy/beam'
		print 'HI systemic velocity = '+str(round(velocity,3))+' km/s'

		return freq, velocity

	def tau_abs(self,scont,sabs):
		'''
		Estimates the optical depth of an absorption line
		INPUT
			scont: continuum flux (float)
			sabs: flux of the absorbed component (negative float)
		OUTPUT
			hi.tau: optical depth of the line	
		'''

		tau=np.log(1-(sabs/scont))
		
		print 'Optical depth = '+str(round(tau,3))

		return tau
		
	def nhi_abs(self,tau,dv):
		'''
		Estimates the column density of the absorption line
		INPUT
			tau: optical depth of the line (float)
			dv: width of the line in km/s
		OUTPUT
			hi.nhi_abs: column density of the absorption line in cm-2	
		'''
		
		nhiabs=self.knhi*self.T*tau*dv
		
		print 'N(HI) = '+str(round(nhiabs,3))+' cm-2'

		return nhiabs

	def beam_area(self, bx,by,z):
		'''
		Estimates the area of the beam of the observations
		INPUT
			bx: x axis of the beam in arcsec (float)
			by: y axis of the beam in arcsec (float)
		OUTPUT
			hi.beamarea: area of the beam in cm2	
		'''

		bxcm=c.ang2lin(bx,z)*1e6
		bycm=c.ang2lin(by,z)*1e6
		beamarea=(bxcm*bycm)*(PC**2)

		print 'Beam Area = '+str(round(beamarea,3))+' cm2'

		return beamarea

	def mhi_abs(self,nhi_abs,area):
		'''
		Estimates the column density of the absorption line
		INPUT
			nhi_abs: column density of the absorption line in cm-2	
			area: area over which the column density is integrated in cm	
			      output of hi.beam_area
		OUTPUT
			hi.mhi_abs: hi mass inferred by the absorption line in Msun
		'''
		
		mhiabs=area*self.mp*nhi_abs/self.msun
		
		print 'M(HI) = '+str(round(mhiabs,3))+' cm-2'	
		
		return mhiabs

				
	def nhi_em(self,s,bx,by,dv):
		'''
		Estimates the column density of the absorption line
		INPUT
			s: flux of the emission line in mJy
			bx: x axis of the beam in arcsec (float)
			by: y axis of the beam in arcsec (float)
			dv: width of the line in km/s
		OUTPUT
			hi.nhi_abs: column density of the absorption line in cm-2	
		'''
		
		bxa=bx/60
		bya=by/60   #convert beam in arcmin
		nhiem=self.nhiem*dv*s/(bxa*bya)
		
		print 'N(HI) = '+str(round(nhi_em,3))+' cm-2'

		return nhi_em

	def mhi_em(self, nhi, area):
		'''
		Estimates the column density of the absorption line
		INPUT
			nhi: column density of the emission line
			area:area over which the column density is integrated in cm	
			      output of hi.beam_area
		OUTPUT
			hi.mhiem: mass inferred from the HI absoprtion line in Msun	
		'''
		
		mhiem=nhi*self.mhi*area/MSUN
		
		print 'M(HI) = '+str(round(mhiem,3))+' cm-2'

		return mhiem

	def mhi_flux(self,z,s,bx,by,pix):
		'''
		Estimates the HI mass from the observed flux, the distance of the
		source, and the resolution of the observations
		(beam and pixel size in arcsec)
		INPUT
			z: redshift of the source (float)
			bx: x axis of the beam in arcsec (float)
			by: y axis of the beam in arcsec (float)
			pix: pixel size in arcsec (float)
		OUTPUT
			hi.mhflux: mass inferred from the flux of the HI emission line	
		'''
		
		dl=c.lum_dist(z)/3.085678e24
		beamcorr=pix**2/(bx*by)
		bla=s*beamcorr
		mhi=self.chi*(dl**2)*bla
		
		return mhi


		
