#!/usr/bin/env python
import os, sys
import yaml

from lmfit import Model
from lmfit.models import GaussianModel
from lmfit.model import save_modelresult

from astropy.io import ascii, fits
from astropy.table import Table, Column
import numpy as np
import numpy.ma as ma

import shutil


import tPlay

tP = tPlay.tplay()



class momplay:

    def makeMoments(self,cfg_par):

        workDir = cfg_par['general']['workdir']

        f = fits.open(workDir+cfg_par['general']['dataCubeName'])
        dd = f[0].header

        lineInfo = tP.openLineList()

        for ii in xrange(0,len(lineInfo['ID'])):
            lineName = str(lineInfo['Name'][ii])
            if '[' in lineName:
                lineName = lineName.replace("[", "")
                lineName = lineName.replace("]", "")

            lineName = lineName+str(int(lineInfo['Wave'][ii]))

            self.moments(lineName,dd,cfg_par['general']['outTableName'])

        return


    def moments(self,lineName,header,outTableName):

        modName = cfg_par['gFit']['modName']
        momModDir = momDir+modName+'/'

        if not os.path.exists(momModDir):
            os.mkdir(momModDir)

        if 'CUNIT3' in header:
            del header['CUNIT3']
        if 'CTYPE3' in header:
            del header['CTYPE3']
        if 'CDELT3' in header:
            del header['CDELT3']
        if 'CRVAL3' in header:  
            del header['CRVAL3']
        if 'CRPIX3' in header:
            del header['CRPIX3'] 
        if 'NAXIS3' in header:
            del header['NAXIS3']

        mom0Head = header.copy()
        mom1Head = header.copy()
        mom2Head = header.copy()
        binHead  = header.copy()


        hdul = fits.open(cfg_par['general']['outTableName'])
        lines = hdul['LineRes_'+cfg_par['gFit']['modName']].data

        tabGen = hdul['BinInfo'].data

        mom0G1 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
        mom1G1 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
        mom2G1 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan

        binMap = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
        
        if modName != 'g1':
            mom0Tot = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
            mom0G2 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
            mom1G2 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
            mom2G2 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
            if modName == 'g3':
                mom0G3 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
                mom1G3 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
                mom2G3 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan


        for i in xrange(0,len(lines['BIN_ID'])):

            match_bin = np.where(tabGen['BIN_ID']==lines['BIN_ID'][i])[0]

            for index in match_bin:
                mom0G1[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['g1_Amp_'+lineName][i]
                mom1G1[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['g1_Centre_'+lineName][i]
                mom2G1[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['g1_sigma_'+lineName][i]
                binMap[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['BIN_ID'][i]

                if modName != 'g1':
                    mom0G2[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['g2_Amp_'+lineName][i]
                    mom1G2[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['g2_Centre_'+lineName][i]
                    mom2G2[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['g2_sigma_'+lineName][i]
                    
                    if modName == 'g3':
                        mom0G3[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['g3_Amp_'+lineName][i]
                        mom1G3[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['g3_Centre_'+lineName][i]
                        mom2G3[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lines['g3_sigma_'+lineName][i]

        binHead['SPECSYS'] = 'topocent'
        binHead['BUNIT'] = 'Flux'
        fits.writeto(momModDir+'binMap.fits',binMap, binHead,overwrite=True)


        mom0Head['SPECSYS'] = 'topocent'
        mom0Head['BUNIT'] = 'Jy/beam.km/s'
        fits.writeto(momModDir+'mom0_g1-'+lineName+'.fits',mom0G1,mom0Head,overwrite=True)

        mom1Head['SPECSYS'] = 'topocent'
        mom1Head['BUNIT'] = 'km/s'
        fits.writeto(momModDir+'mom1_g1-'+lineName+'.fits',mom1G1,mom1Head,overwrite=True)

        mom2Head['SPECSYS'] = 'topocent'
        mom2Head['BUNIT'] = 'km/s'
        fits.writeto(momModDir+'mom2_g1-'+lineName+'.fits',mom2G1,mom2Head,overwrite=True)
        
        if modName != 'g1':
            fits.writeto(momModDir+'mom0_g2-'+lineName+'.fits',mom0G2,mom0Head,overwrite=True)
            fits.writeto(momModDir+'mom1_g2-'+lineName+'.fits',mom1G2,mom1Head,overwrite=True)
            if modName == 'g2':
                fits.writeto(momModDir+'mom0_tot-'+lineName+'.fits',mom0G1+mom0G2,mom0Head,overwrite=True)
            if modName == 'g3':
                fits.writeto(momModDir+'mom0_g3-'+lineName+'.fits',mom0G3,mom0Head,overwrite=True)
                fits.writeto(momModDir+'mom1_g3-'+lineName+'.fits',mom1G3,mom1Head,overwrite=True)
                fits.writeto(momModDir+'mom2_g3-'+lineName+'.fits',mom2G3,mom2Head,overwrite=True)
                fits.writeto(momModDir+'mom0_tot-'+lineName+'.fits',mom0G1+mom0G2+mom0G3,mom0Head,overwrite=True)

        return

    def makeLineRatioMaps(self,cfg_par):

        workDir = cfg_par['general']['workdir']

        f = fits.open(workDir+cfg_par['general']['dataCubeName'])
        dd = f[0].header

        #lineInfo = self.openLineList()

        #for ii in xrange(0,len(lineInfo['ID'])):
        #    lineName = str(lineInfo['Name'][ii])
        #    if '[' in lineName:
        #        lineName = lineName.replace("[", "")
        #        lineName = lineName.replace("]", "")

            #lineName = lineName+str(int(lineInfo['Wave'][ii]))

        self.momLineRatio(dd,cfg_par['general']['outTableName'])

        return


    def momLineRatio(self,header,outTableName):


        modName = cfg_par['gFit']['modName']
        momModDir = cfg_par['general']['momDir']+modName+'/'

        if not os.path.exists(momModDir):
            os.mkdir(momModDir)

        if 'CUNIT3' in header:
            del header['CUNIT3']
        if 'CTYPE3' in header:
            del header['CTYPE3']
        if 'CDELT3' in header:
            del header['CDELT3']
        if 'CRVAL3' in header:  
            del header['CRVAL3']
        if 'CRPIX3' in header:
            del header['CRPIX3'] 
        if 'NAXIS3' in header:
            del header['NAXIS3']

        lineMapHead = header.copy()
        
        hdul = fits.open(cfg_par['general']['outTableName'])
        lineBPT = hdul['BPT_'+cfg_par['gFit']['modName']].data

        tabGen = hdul['BinInfo'].data

        numCols = len(lineBPT.dtype.names)
        if modName == 'g2':
            numCols = (numCols-1)/3
        if modName == 'g3':
            numCols = (numCols-1)/4
        print numCols

        for i in xrange(1,numCols+1):
            lineMapG1 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
        
            if modName != 'g1':
                lineMapToT = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
                lineMapG2 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan
            
            if modName == 'g3':
            
                lineMapG3 = np.zeros([header['NAXIS2'],header['NAXIS1']])*np.nan

            for j in xrange(0,len(lineBPT['BIN_ID'])):

                match_bin = np.where(tabGen['BIN_ID']==lineBPT['BIN_ID'][j])[0]
                for index in match_bin:
                    lineMapG1[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lineBPT[j][i]

                    if modName != 'g1':                        
                        lineMapToT[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lineBPT[j][i+numCols*2]
                        lineMapG2[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lineBPT[j][i+numCols]

                        if modName == 'g3':
                            lineMapG3[int(tabGen['PixY'][index]),int(tabGen['PixX'][index])] = lineBPT[j][i+numCols+2]



            lineMapHead['BUNIT'] = 'Flux'
            fits.writeto(momModDir+'BPT-'+str(lineBPT.dtype.names[i])+'.fits',lineMapG1,lineMapHead,overwrite=True)
            
            if modName != 'g1':
                fits.writeto(momModDir+'BPT-'+str(lineBPT.dtype.names[i+numCols])+'.fits',lineMapG2,lineMapHead,overwrite=True)
                fits.writeto(momModDir+'BPT-'+str(lineBPT.dtype.names[i+numCols*2])+'.fits',lineMapToT,lineMapHead,overwrite=True)

                if modName == 'g3':
                    fits.writeto(momModDir+'BPT-'+str(lineBPT.dtype.names[i+numCols+2])+'.fits',lineMapG3,lineMapHead,overwrite=True)

        return