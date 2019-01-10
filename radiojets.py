#!/usr/bin/env python

import sys, string, os
import numpy as np

from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u

import casacore.tables as tables
import logging 


class radiojets:

	def __init__(self):

		self.rootdir = os.getcwd()+'/'


        	#set logger
        	logging.basicConfig(level=logging.INFO)
		self.logger = logging.getLogger(__name__)  


	def orientationAngle(self,value,spIndex):

		rpow = np.power(value,1./(2+spIndex))
		

		orAngle = np.divide(rpow-1. , 1+rpow)

		return orAngle

	def brightnessRatio(self,a,b):

		ratio =  np.divide(a,b)

		return ratio	



#-------------------------------#
#             MAIN              #
#-------------------------------#
rjet=radiojets()

ops=sys.argv[1]

if 'help' in sys.argv or '-h' in sys.argv:
    print 'Run as follows'
    print 'radiojets.py -ops operation OR -v values '
    sys.exit()
else: arg=sys.argv

ops=arg[arg.index('-op')+1]

if len(arg) == 6:
	v_a = np.float(arg[arg.index('-v')+1])
	v_b = np.float(arg[arg.index('-v')+2])
else:
	v_a = np.float(arg[arg.index('-v')+1])

if ops == 'bRatio' or ops == 'brightnessRatio' or ops == 'bR':

	ratio = rjet.brightnessRatio(v_a,v_b)
	print ratio

elif ops == 'orAngle' or ops == 'orientationAngle' or ops == 'oA':

	orAngle = rjet.orientationAngle(v_a,v_b)
	print orAngle
	print orAngle/np.pi*180.



