# Import modules
import sys, string, os
import numpy as np
import casacore.tables as tables
import logging 


class msinfo:

	def __init__(self):

		self.rootdir = os.getcwd()+'/'
    		self.filename=sys.argv[1]

        	#set logger
        	logging.basicConfig(level=logging.INFO)
		self.logger = logging.getLogger(__name__)  
		
    	def time_chunk(self):

		self.logger.info("\t ...  Observing time Info ... \n")

		t=tables.table(self.msfile)
		self.time = t.getcol('TIME')
		t.close()

		starttime= self.time[0]
		endtime=self.time[-1]
		time_chunk = float(cfg_par['rfi']['chunks']['time_step'])*60.

		times=np.arange(starttime,endtime+time_chunk*1.,time_chunk)

		startdate=Time(starttime/3600./24.,format='mjd',scale='utc')
		cfg_par['rfi']['startdate'] = startdate
		startdate.format='iso' 
		startdate.subformat='date_hm'       

		enddate=Time(endtime/3600./24.,format='mjd',scale='utc')
		cfg_par['rfi']['enddate'] = enddate

		enddate.format='iso'        
		enddate.subformat='date_hm'       

		self.logger.info('\t Start date: {0:%y}{0:%b}{0:%d}:{0:%X}'.format(startdate.datetime))
		self.logger.info('\t End date  : {0:%y}{0:%b}{0:%d}:{0:%X} \n\n'.format(enddate.datetime))
        
		return times,startdate,enddate
  	
	def load_from_ms(self):
		'''
		Loads important columns from MS file
		From MS: 
		Field_ID, Data, Flag, Antenna1, Antenna2,
		From MS/ANTENNA:
		Position,Name
		From MS/SPECTRAL_WINDOW
		Chan_width, Chan_freq
		'''

		self.logger.info("\t ... Fields, Antennas & Bandwidth Info ...\n")

		self.msfile = self.filename
		
		fields=tables.table(self.msfile+'/FIELD')
		self.fieldNames = fields.getcol('NAME')
		
		for i inxrange(0,len(self.fieldNames)):
			self.coords=fields.getcol('REFERENCE_DIR')
			self.coords =self.coords*180./np.pi
			ra,dec = SkyCoord(self.coords[self.fieldNames,:,0]*u.degree, self.coords[self.fieldNames,:,1]*u.degree,  unit=(u.deg, u.deg))

			self.logger.info("\tField with name {0:s} (Field ID = {1:d}): ".format(self.FieldName[i],i,ra,dec))
		    #self.logger.info("\tCoordinates {}".format(selectFieldName,self.selectFieldID))

		antennas = tables.table(self.msfile +'/ANTENNA')
		self.ant_pos = np.array(antennas.getcol('POSITION'))
		self.ant_wsrtnames = np.array(antennas.getcol('NAME'))

		self.ant_names = np.arange(0,self.ant_wsrtnames.shape[0],1)
		self.nant = len(self.ant_names)

		#logging
		self.logger.info("\tTotal number of antennas:\t"+str(self.nant))
		self.logger.info("\tAntenna names:\t\t"+str(self.ant_names))

		antennas.close()

		spw=tables.table(self.msfile+'/SPECTRAL_WINDOW')
		self.channelWidths=spw.getcol('CHAN_WIDTH')

		self.channelFreqs=spw.getcol('CHAN_FREQ')
		self.chan_widths = self.channelWidths[0][0]
		self.lowfreq = float(self.channelFreqs[0][0])
		self.highfreq = float(self.channelFreqs[-1][-1])

		spw.close()

		self.logger.info("\tChannel Width [kHz]:\t"+str(cfg_par['rfi']['chan_widths']/1e3))
		self.logger.info("\tStart         [GHz]:\t"+str(cfg_par['rfi']['lowfreq']/1e9))
		self.logger.info("\tEnd           [GHz]:\t"+str(cfg_par['rfi']['highfreq']/1e9)+'\n')


		#determine start and end date
		times_tm, start_tmp, end_tmp = self.time_chunk(cfg_par)


		self.logger.info("\t ... info from MS file loaded  \n\n")

		return 0
