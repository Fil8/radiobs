# Import modules
import sys, string, os
import numpy as np
import casacore.tables as tables
import logging 


class msinfo:

	def __init__(self):

	  self.rootdir = os.getcwd()+'/'
    self.filename=sys.argv[1]
  
  def load_from_ms(self,cfg_par,times=0):
    '''
    Loads important columns from MS file
    From MS: 
    Field_ID, Data, Flag, Antenna1, Antenna2,
    From MS/ANTENNA:
    Position,Name
    From MS/SPECTRAL_WINDOW
    Chan_width, Chan_freq
    '''

    self.logger.info("\t ... Field, Antenna & Bandwidth Info ...\n")

    self.msfile = cfg_par['general']['msfullpath']
    self.aperfi_badant = cfg_par['rfi']['bad_antenna'] 
    self.selectFieldID = cfg_par['general']['field']



    fields=tables.table(self.msfile+'/FIELD')
    self.fieldNames = fields.getcol('NAME')
    selectFieldName= self.fieldNames[int(self.selectFieldID)]
    cfg_par['general']['fieldname'] = selectFieldName
    self.coords=fields.getcol('REFERENCE_DIR')
    self.coords =self.coords*180./np.pi
    cfg_par['rfi']['coords'] = SkyCoord(self.coords[self.selectFieldID,:,0]*u.degree, self.coords[self.selectFieldID,:,1]*u.degree,  unit=(u.deg, u.deg))

    self.logger.info("\tField with name {0:s} (Field ID = {1:d})".format(selectFieldName,self.selectFieldID))
    #self.logger.info("\tCoordinates {}".format(selectFieldName,self.selectFieldID))

    antennas = tables.table(self.msfile +'/ANTENNA')
    self.ant_pos = np.array(antennas.getcol('POSITION'))
    self.ant_wsrtnames = np.array(antennas.getcol('NAME'))

    self.ant_names = np.arange(0,self.ant_wsrtnames.shape[0],1)
    self.nant = len(self.ant_names)

    #logging
    cfg_par['rfi']['nant'] = self.nant
    cfg_par['rfi']['ant_names'] = self.ant_names

    self.logger.info("\tTotal number of antennas:\t"+str(self.nant))
    self.logger.info("\tAntenna names:\t\t"+str(self.ant_names))

    antennas.close()

    spw=tables.table(self.msfile+'/SPECTRAL_WINDOW')
    self.channelWidths=spw.getcol('CHAN_WIDTH')

    self.channelFreqs=spw.getcol('CHAN_FREQ')
    cfg_par['rfi']['chan_widths'] = self.channelWidths[0][0]
    cfg_par['rfi']['lowfreq'] = float(self.channelFreqs[0][0])
    cfg_par['rfi']['highfreq'] = float(self.channelFreqs[-1][-1])

    spw.close()

    self.logger.info("\tChannel Width [kHz]:\t"+str(cfg_par['rfi']['chan_widths']/1e3))
    self.logger.info("\tStart         [GHz]:\t"+str(cfg_par['rfi']['lowfreq']/1e9))
    self.logger.info("\tEnd           [GHz]:\t"+str(cfg_par['rfi']['highfreq']/1e9)+'\n')


    #determine start and end date
    times_tm, start_tmp, end_tmp = rfiST.time_chunk(cfg_par)

    t=tables.table(self.msfile)

    if times !=0:
        value_end = times[1]
        value_start = times[0]

        altaz = rfiST.alt_az(cfg_par,times[0])
        cfg_par['rfi']['altaz'] = altaz

        t2 = tables.taql('select from $t where TIME < $value_end and TIME>$value_start')
        self.fieldIDs=t2.getcol('FIELD_ID')
        self.ant1 = t2.getcol('ANTENNA1')
        self.ant2 = t2.getcol('ANTENNA2')

        selection  = self.fieldIDs==self.selectFieldID
        selection *= self.ant1!=self.ant2
        if np.sum(selection) !=0. :
            if cfg_par['rfi']['RFInder_mode'] == 'rms_clip':
                self.vis  = t2.getcol('DATA')[selection]
            self.flag = t2.getcol('FLAG')[selection]
            self.interval = t2.getcol('INTERVAL')[selection]
            empty_table=0     
        else:
            self.logger.warning('\t ### Table of selected interval is empty ')
            self.logger.warning('\t     Correct noise_measure_edges in rfi of parameter file ###')
            empty_table=1
        t2.close()

    else:

        self.fieldIDs=t.getcol('FIELD_ID')
        self.ant1=t.getcol('ANTENNA1')
        self.ant2=t.getcol('ANTENNA2')

        #select from cross correlations of the correct field
        selection=self.fieldIDs==self.selectFieldID
        selection*=self.ant1!=self.ant2

        altaz = rfiST.alt_az(cfg_par,times_tm[0])
        cfg_par['rfi']['altaz'] = altaz

        if np.sum(selection) !=0. :
            if cfg_par['rfi']['RFInder_mode'] == 'rms_clip':
                self.vis  = t.getcol('DATA')[selection]
            self.flag = t.getcol('FLAG')[selection]
            self.interval = t.getcol('INTERVAL')[selection]
            empty_table=0     
        else:
            self.logger.warning('\t ### Table of selected interval is empty ')
            self.logger.warning('\t     Correct noise_measure_edges in rfi of parameter file ###')
            empty_table=1


    t.close()

    if not self.aperfi_badant:
        nrbadant =len(int(self.aperfi_badant))
    else:
        nrbadant = 0.

    #visibilities over all times per baseline
    #estimate noise
    if empty_table!=1:
        #determine number of baselines
        nrAnt=np.unique(np.concatenate((self.ant1,self.ant2))).shape[0]
        nrBaseline=(nrAnt-nrbadant)*(nrAnt-nrbadant-1)/2        
        cfg_par['rfi']['number_baseline'] = nrBaseline
        rfiST.predict_noise(cfg_par,self.channelWidths,self.interval,self.flag)
        cfg_par['rfi']['vis_alltimes_baseline'] = self.flag.shape[0]/nrBaseline

    self.logger.info("\t ... info from MS file loaded  \n\n")

    return empty_table
