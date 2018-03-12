__author__ = "Filippo Maccagni"
__copyright__ = "Fil8"
__email__ = "filippo.maccagni@gmail.com"

import numpy as np

class kk:
	'''
	Useful constants called by other radiobs classes as kk.constant
	'''

	#converstions
	RAD2DEG=180./np.pi
	PC=3.08567758E18    #cm
	JANSKY=1e-23        #erg/scm2Hz
	
	#constants
	HI=1.42040575177e+09 #Hz
	TSPIN=100            #K
	C=2.99792458E10     #cm/s
	G=6.6742E-08        #cm3kg-1s-1 
	MSUN=1.98855e33      #g
	MHI=1.6749E-24       #g
	CHI=2.36E5
	MP=1.67492728E-24   #g
	SIGMAT=6.66524E-25  #cm2
	
	#convert to column density of HI
	KnhiABS = 1.8216E18
	KnhiEM = 3.1E17
	
	#magnitudes & factors for Tully-Fisher
	M26 = -23.33
	STULLY = -9.64

	#cosmological constants
	h0 = 70
	omega_l = 0.7 
	omega_m = 0.3