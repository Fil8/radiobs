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
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib import gridspec
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, LogLocator
from matplotlib import transforms as mtransforms
from matplotlib.ticker import LogFormatter 
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

class gplay:
    

    def __init__(self,file=None):

        #self.rootdir = os.getcwd()+'/'
        self.C = 2.99792458e8

        # get directories
        GFIT_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        GFIT_DIR = GFIT_PATH+'/ants/'
        sys.path.append(os.path.join(GFIT_PATH, 'radiobs'))
        
        file_default = GFIT_DIR + 'gPlay_default.yaml'

        if file != None:
            cfg = open(file)
        else:
            cfg = open(file_default)

        
        self.cfg_par = yaml.load(cfg)
   
        runDir = self.cfg_par['general']['workdir']+self.cfg_par['general']['runName']
        if not os.path.exists(runDir):
            os.mkdir(runDir)
        self.cfg_par['general']['runNameDir'] = runDir
        #if self.cfg_par['general']['verbose'] == True:
        #    print yaml.dump(self.cfg_par)


        cfg.close()

        return

    def lambdaVRad(self,lambdaWave,lambdaRest):

        # wavelenght in Angstrom
        vRad = ((lambdaWave-lambdaRest)/lambdaRest)*self.C/1e3

        #velocity in km/s
        return vRad

    def vRadLambda(self,vRad,lambdaRest):
        
        # velocity in kms
        lambdaWave = (vRad*1e3*lambdaRest*1e-10)/(self.C)/1e-10 +lambdaRest

        # wavelenght in Angstrom
        return lambdaWave

    def modDef(self,modName):

        
        Gmod = GaussianModel()

        gauss1 = GaussianModel(prefix='g1_ln'+str(i)+'_')

        pars = gauss1.make_params()

        pars['g1_ln'+0+'_'+'center'].set(value=np.argmax[y[indexMin:indexMax]],
            min=waveMin,max=waveMax)
        pars['g1_ln'+str(i)+'_'+'height'].set(value=ampIn1)
        pars['g1_ln'+str(i)+'_'+'sigma'].set(expr='g1_sigma')


        gauss2 = GaussianModel(prefix='g2_')
        pars.update(gauss2.make_params())

        pars['g2_center'].set(value=-10, min=wMin, max=wMax)
        pars['g2_sigma'].set(value=100, min=0)
        pars['g2_amplitude'].set(value=2000, min=1)

        cenIn3 = cenIn1 + dltV13
        sigmaIn3 = pars['g1_ln'+str(i)+'_'+'sigma'] + dltSgm12
        ampIn3 = ampIn1 + dltAmp13

        gauss3 = GaussianModel(prefix='g3_')
        pars.update(gauss3.make_params())

        pars['g3_center'].set(value=-20, min=wMin, max=wMax)
        pars['g3_sigma'].set(value=200, min=0)
        pars['g3_amplitude'].set(value=2000, min=1)

        gName1 = ['g1']
        gName2 = ['g1', 'g2']   
        gName3 = ['g1', 'g2', 'g3']
        
        mod1 = gauss1 
        mod2 = gauss1+gauss2
        mod3 = gauss1+gauss2+gauss3

        if modName == 'g1':
            return mod1,gName1
        if modName == 'g2':
            return mod2,gName2
        if modName == 'g3':
            return mod3,gName3

    def lineModDef(self,wave,y,lineInfo):

        gName = self.cfg_par['gFit']['modName']
        for i in xrange(0,len(lineInfo['ID'])):
            
            #waveMin =  np.log(lineInfo['Wave'][i] - lineInfo['lineRangeAng'][i])
            #waveMax =  np.log(lineInfo['Wave'][i] + lineInfo['lineRangeAng'][i])
           
            waveAmpIn1Min = np.log(lineInfo['Wave'][i]-lineInfo['cenRangeAng'][i])
            indexMin = int(np.where(abs(wave-waveAmpIn1Min)==abs(wave-waveAmpIn1Min).min())[0]) 
    
            waveAmpIn1Max = np.log(lineInfo['Wave'][i]+lineInfo['cenRangeAng'][i])
            indexMax = int(np.where(abs(wave-waveAmpIn1Max)==abs(wave-waveAmpIn1Max).min())[0])
                       
            Gmod = GaussianModel()

            gauss1 = GaussianModel(prefix='g1ln'+str(i)+'_')

            sigmaIn1 = lineInfo['deltaSigmaAng_In1'][i]
            ampIn1 = np.max(y[indexMin:indexMax])*max(2.220446049250313e-16, sigmaIn1)/0.3989423
            smallWave = wave[indexMin:indexMax]
            cenIn1 = smallWave[np.argmax(y[indexMin:indexMax])]

            if i == 0:
                pars = gauss1.make_params()
                pars['g1ln'+str(i)+'_'+'sigma'].set(value=sigmaIn1,
                    min=sigmaIn1/5.,max=sigmaIn1*5.)
                pars['g1ln'+str(i)+'_'+'center'].set(value=cenIn1,
                min=waveAmpIn1Min,max=waveAmpIn1Max)
                pars['g1ln'+str(i)+'_'+'amplitude'].set(value=ampIn1,min=0,max=None)
                mod = gauss1
          
            else:
            
                pars.update(gauss1.make_params())    
                pars['g1ln'+str(i)+'_'+'center'].set(value=cenIn1,
                        min=waveAmpIn1Min,max=waveAmpIn1Max)
                pars['g1ln'+str(i)+'_'+'amplitude'].set(value=ampIn1,min=0,max=None)
                if self.cfg_par['gFit']['fixSigma'] == True:
                    pars['g1ln'+str(i)+'_'+'sigma'].set(expr='g1ln'+str(0)+'_'+'sigma')
                else:
                    pars['g1ln'+str(i)+'_'+'sigma'].set(value=sigmaIn1,
                    min=sigmaIn1/5.,max=sigmaIn1*5.)    

                mod += gauss1            
            
            if gName != 'g1':

                Gmod = GaussianModel()
                gauss2 = GaussianModel(prefix='g2ln'+str(i)+'_')
                pars.update(gauss2.make_params())

                if i == 0:
                    sigmaIn2 = pars['g1ln'+str(i)+'_'+'sigma'] +lineInfo['deltaSigmaAng_12'][i]

                    pars['g2ln'+str(i)+'_'+'sigma'].set(value=sigmaIn2,min=sigmaIn2/5.,max=sigmaIn2*5.)
                elif self.cfg_par['gFit']['fixSigma'] == True:
                    pars['g2ln'+str(i)+'_'+'sigma'].set(expr='g2ln'+str(0)+'_'+'sigma')
                else:
                    pars['g2ln'+str(i)+'_'+'sigma'].set(value=sigmaIn2,min=sigmaIn2/5.,max=sigmaIn2*5.)
                
                ampIn2 = ampIn1*self.cfg_par['gFit']['dltAmp12']               
                cenIn2Pos = cenIn1 + lineInfo['deltaVAng_12'][i]
                cenIn2Neg = cenIn1 - lineInfo['deltaVAng_12'][i]
              
                if self.cfg_par['gFit']['fixCentre'] == True and i >0:

                    pars['g2ln'+str(i)+'_'+'center'].set(value=cenIn2Pos,
                        min=waveAmpIn1Min-lineInfo['deltaVAng_12'][i],max=waveAmpIn1Max+lineInfo['deltaVAng_12'][i],vary=True)
                    
                    pars.add(name='g2ln'+str(i)+'Split_'+'center', expr='g2ln'+str(0)+'_'+'center - g1ln'+str(0)+'_'+'center')
                    pars.add(name='g2ln'+str(i)+'Pos_'+'center', value=cenIn2Pos,
                        max=waveAmpIn1Max+lineInfo['deltaVAng_12'][i],min=cenIn1, vary=True)
                    pars.add(name='g2ln'+str(i)+'Neg_'+'center', value=cenIn2Neg,max=cenIn1,
                        min=waveAmpIn1Min-lineInfo['deltaVAng_12'][i], vary=True)

                    pars['g2ln'+str(i)+'_'+'center'].set(expr='g2ln'+str(i)+'Pos_center if g2ln'+str(i)+'Split_center >= 0 else g2ln'+str(i)+'Neg_center' )                    

                else:
                    pars['g2ln'+str(i)+'_'+'center'].set(value=cenIn2Pos,
                        min=waveAmpIn1Min-lineInfo['deltaVAng_12'][i],max=waveAmpIn1Max+lineInfo['deltaVAng_12'][i])                   

                pars['g2ln'+str(i)+'_'+'amplitude'].set(value=ampIn2,min=0,max=None)
                

                mod += gauss2


                if gName == 'g3':

                    Gmod = GaussianModel()

                    gauss3 = GaussianModel(prefix='g3ln'+str(i)+'_')

                    pars.update(gauss3.make_params())

                    
                    if i == 0:
                        sigmaIn3 = pars['g1ln'+str(i)+'_'+'sigma'] + lineInfo['deltaSigmaAng_13'][i]
                        pars['g3ln'+str(i)+'_'+'sigma'].set(value=sigmaIn3,min=sigmaIn3/5.,max=sigmaIn3*5.)
                    elif self.cfg_par['gFit']['fixSigma'] == True:
                        pars['g3ln'+str(i)+'_'+'sigma'].set(expr='g3ln'+str(0)+'_'+'sigma')
                    else:
                        pars['g3ln'+str(i)+'_'+'sigma'].set(value=sigmaIn3,min=sigmaIn3/5.,max=sigmaIn3*5.)

                    ampIn3 = ampIn1*self.cfg_par['gFit']['dltAmp13']
                    
                    cenIn3Pos = cenIn1 + lineInfo['deltaVAng_13'][i]
                    cenIn3Neg = cenIn1 - lineInfo['deltaVAng_13'][i]
                    
                    pars['g3ln'+str(i)+'_'+'center'].set(value=cenIn3Pos,
                        min=waveAmpIn1Min-lineInfo['deltaVAng_13'][i],max=waveAmpIn1Max+lineInfo['deltaVAng_13'][i])
                
                    if self.cfg_par['gFit']['fixCentre'] == True and i >0:

                        pars['g3ln'+str(i)+'_'+'center'].set(value=cenIn3Pos,
                            min=waveAmpIn1Min-lineInfo['deltaVAng_13'][i],max=waveAmpIn1Max+lineInfo['deltaVAng_13'][i],vary=True)
                        
                        pars.add(name='g3ln'+str(i)+'Split_'+'center', expr='g3ln'+str(0)+'_'+'center - g1ln'+str(0)+'_'+'center')
                        pars.add(name='g3ln'+str(i)+'Pos_'+'center', value=cenIn3Pos,
                            max=waveAmpIn1Max+lineInfo['deltaVAng_13'][i],min=cenIn1, vary=True)
                        pars.add(name='g3ln'+str(i)+'Neg_'+'center', value=cenIn3Neg,max=cenIn1,
                            min=waveAmpIn1Min-lineInfo['deltaVAng_13'][i], vary=True)

                        pars['g3ln'+str(i)+'_'+'center'].set(expr='g3ln'+str(i)+'Pos_center if g3ln'+str(i)+'Split_center >= 0 else g3ln'+str(i)+'Neg_center' )                    

                    else:
                        pars['g3ln'+str(i)+'_'+'center'].set(value=cenIn3Pos,
                            min=waveAmpIn1Min-lineInfo['deltaVAng_13'][i],max=waveAmpIn1Max+lineInfo['deltaVAng_13'][i])                   

                    pars['g3ln'+str(i)+'_'+'amplitude'].set(value=ampIn3,min=0,max=None)

                    mod += gauss3

        #pars.pretty_print()
        return mod,pars



    def gFit(self):
        key = 'general'
        workDir = self.cfg_par[key]['workdir']
        #open line lineList
        lineInfo = self.openLineList(workDir+self.cfg_par[key]['lineListName'])

        #open table for bins
        wave,xAxis,yAxis,pxSize,noiseBin, vorBinInfo = self.openTablesPPXF(workDir+self.cfg_par[key]['tableBinName'],
            workDir+self.cfg_par[key]['tableSpecName'])
        #open datacube
        f = fits.open(workDir+self.cfg_par[key]['dataCubeName'])
        hh = f[0].header
        dd = f[0].data
        

        runNameDir = workDir+self.cfg_par[key]['runName']
        self.cfg_par[key]['runNameDir'] = runNameDir
        if not os.path.exists(runNameDir):
            os.mkdir(runNameDir)
        #call model

        #define x-axis array
        #wave = ((np.linspace(1, dd.shape[0], dd.shape[0]) - hh['CRPIX3']) * hh['CDELT3'] + hh['CRVAL3'])
        lambdaMin = np.log(self.cfg_par['gFit']['lambdaMin'])
        lambdaMax = np.log(self.cfg_par['gFit']['lambdaMax'])


        idxMin = int(np.where(abs(wave-lambdaMin)==abs(wave-lambdaMin).min())[0]) 
        idxMax = int(np.where(abs(wave-lambdaMax)==abs(wave-lambdaMax).min())[0] )# If you want the index of the element of array (array) nearest to the the given number (num)

        # ?
        diffusion = 1e-3
        binID = np.zeros([dd.shape[1],dd.shape[2]])*np.nan

        #loop over X,Y of datacube -> why not over single Voronoi bins?
        for j in xrange(205,206):
            for i in xrange(250,251):
        #for j in xrange(0,dd.shape[1]):
        #    for i in xrange(0,dd.shape[2]):
                
                y = dd[idxMin:idxMax,j,i]

                waveCut = wave[idxMin:idxMax]

                #check if spectrum is not empty                   
                if np.sum(y)>0:


                    gMod,gPars = self.lineModDef(waveCut,y,lineInfo)


                    # identify voronoi bin
                    xVal = xAxis[i]
                    yVal = yAxis[j]

                    index = np.where((vorBinInfo['X'] < (xVal+pxSize/2.)) & 
                    ((xVal-pxSize/2.) < vorBinInfo['X']) & (vorBinInfo['Y'] < (yVal+pxSize/2.)) & 
                    ((yVal-pxSize/2.) < vorBinInfo['Y']))

                    binIDName = vorBinInfo['BIN_ID'][index]
                    IDName = vorBinInfo['ID'][index]
                    xx = vorBinInfo['X'][index]
                    yy = vorBinInfo['Y'][index]
                    PixX = i
                    PixY = j
  
                    binFitName =runNameDir+str(binIDName)+'fit.sav'

                    singleVorBinInfo = {'ID':IDName,'BIN_ID':binIDName, 'X':xx, 'Y':yy, 'PixX':PixX, 'PixY': PixY}
 
                    #check if it is first time in bin
                    if binIDName not in binID[:,:]:
                        binID[j,i] = binIDName
                        noiseVec = noiseBin[binIDName][:] [0] 

                        # FIT
                        result = gMod.fit(y, gPars, x=waveCut)
                        #print(result.fit_report())

                        
                        #plot Fit
                        #self.plotSpecFit(waveCut, y,result,noiseVec[idxMin:idxMax],i,j,lineInfo)
                        #self.plotLineZoom(waveCut, y,result,noiseVec[idxMin:idxMax],i,j,lineInfo)

                        #save Fit for fast upload 

                        #save_modelresult(result,binFitName)

                        #mom0[j,i] = fitRes['amplitude']
                        #mom1[j,i] = fitRes['center']

                    #else:

#                        result = load_modelresult(binFitName)
#                        self.loadFitFromRightBin(binNumber)



                    #self.updateFitTable(result,gName)
                        #self.writeOutputTable(lineInfo,singleVorBinInfo,result) 
                                                
                else:
                    pass

        print('''\t+---------+\n\t fits done\n\t+---------+''')
    
        return lineInfo,singleVorBinInfo,result
    
    def openLineList(self,lineList):
        
        lineInfo = ascii.read(lineList) 

        #mask line list 
        index = np.where(lineInfo['Fit'] == 0)

        indexLines = np.where(lineInfo['Fit'] == 1)
        #idxMin = 

        fltr =  np.array(index)[0]
        lineInfo.remove_rows(list(fltr))

        lenTable = len(lineInfo['ID'])
        dltSigmaIn1Ang = np.zeros([lenTable])

        dltV12Ang = np.zeros([lenTable])
        dltSigma12Ang = np.zeros([lenTable])

        dltV13Ang = np.zeros([lenTable])
        dltSigma13Ang = np.zeros([lenTable])
        lineRange = np.zeros([lenTable])
        cenRange = np.zeros([lenTable])

        for i in xrange(0,lenTable):

            lambdaRest = lineInfo['Wave'][i]
            
            lineRange[i] = self.vRadLambda(lineInfo['lineRange'][i],
                lambdaRest)-lambdaRest    

            cenRange[i] = self.vRadLambda(lineInfo['cenRange'][i],
                lambdaRest)-lambdaRest    

            deltaV12 = np.log(self.vRadLambda(self.cfg_par['gFit']['dltV12'],
                lambdaRest))
            deltaV12 -= np.log(lambdaRest)
            
            deltaV13 =np.log(self.vRadLambda(self.cfg_par['gFit']['dltV13'],
                lambdaRest))
            deltaV13 -= np.log(lambdaRest)


            deltaSigmaIn1 = np.log(self.vRadLambda(self.cfg_par['gFit']['sigmaIn1'],
                lambdaRest))
            deltaSigmaIn1 -= np.log(lambdaRest)
            
            deltaSigma12 = np.log(self.vRadLambda(self.cfg_par['gFit']['dltSigma12'],
                lambdaRest))
            deltaSigma12 -=  np.log(lambdaRest)
            
            deltaSigma13 = np.log(self.vRadLambda(self.cfg_par['gFit']['dltSigma13'],
                lambdaRest))
            deltaSigma13 -= np.log(lambdaRest)


            dltSigmaIn1Ang[i] = deltaSigmaIn1
            dltV12Ang[i] = deltaV12            
            dltSigma12Ang[i] = deltaSigma12

            dltV13Ang[i] = deltaV13
            dltSigma13Ang[i] = deltaSigma13

        dltSigmaIn1Col = Column(name='deltaSigmaAng_In1', data=dltSigmaIn1Ang)        
        dltV12Col = Column(name='deltaVAng_12', data=dltV12Ang)
        dltSigma12Col = Column(name='deltaSigmaAng_12', data=dltSigma12Ang)
        dltV13Col = Column(name='deltaVAng_13', data=dltV13Ang)
        dltSigma13Col = Column(name='deltaSigmaAng_13', data=dltSigma13Ang)
        lineRangeCol = Column(name='lineRangeAng', data=lineRange)
        cenRangeCol = Column(name='cenRangeAng', data=cenRange)


        lineInfo.add_column(dltSigmaIn1Col)
        lineInfo.add_column(dltV12Col)
        lineInfo.add_column(dltSigma12Col)
        lineInfo.add_column(dltV13Col)
        lineInfo.add_column(dltSigma13Col)
        lineInfo.add_column(lineRangeCol)
        lineInfo.add_column(cenRangeCol)


        return lineInfo

    def openTablesPPXF(self,tableBin,tableSpec):
        crPix1=200 # E 
        crPix2=125
      
        tab = fits.open(tableBin)
        head = tab[0].header
        headTab = tab[1].header
        dataTab = tab[1].data    
        head['CRPIX1'] = crPix1
        head['CRPIX2'] = crPix2 
        
        xMin = np.min(dataTab['X'])
        xMax = np.max(dataTab['X'])

        shapeX = (xMax-xMin)/head['PIXSIZE']

        yMin = np.min(dataTab['Y'])
        yMax = np.max(dataTab['Y'])

        shapeY = (yMax-yMin)/head['PIXSIZE']

        xAxis = (np.linspace(1, shapeX+1, shapeX+1)-head['CRPIX1']+1) *head['PIXSIZE']
        yAxis = (np.linspace(1, shapeY+1, shapeY+1)-head['CRPIX2']+1) *head['PIXSIZE']
        tab = fits.open(tableSpec)
        dataSpec = tab[1].data
        specExp = tab[2].data
        wave = [item for t in specExp for item in t] 

        noiseBin = dataSpec['ESPEC']

        pxSize = head['PIXSIZE']

        return wave,xAxis,yAxis,pxSize,noiseBin,dataTab

    def updateFitTable(self,result,gName):


        for k in xrange (0,len(gName)):

            amp = fitRes[gName[k]+'_amplitude']
            ctr = fitRes[gName[k]+'_center']
            sig = fitRes[gName[k]+'_sigma']

            fwhm = fitRes[gName[k]+'_fwhm']
            height = fitRes[gName[k]+'_height']


            amp_err = result.params[gName[k]+'_amplitude'].stderr
            sig_err = result.params[gName[k]+'_sigma'].stderr
            cen_err = result.params[gName[k]+'_center'].stderr            

            tabErr = {gName[k]+'_amplitude':amp, gName[k]+'_center': ctr, gName[k]+'_sigma': sig,
                      gName[k]+'_fwhm' : fwhm, gName[k]+'_height' : height,
                      gName[k]+'_centreErr':cen_err, gName[k]+'_ampErr': amp_err, 
                      gName[k]+'_sigmaErr': sig_err}

            tabInfo.update(tabErr)

        aic = result.aic
        bic = result.bic
        redchi = result.redchi
        success = result.success
        ndata = result.ndata
        nvarys = result.nvarys
        nfev = result.nfev

        tabFit = {'aic':aic,'bic': bic, 'ndata': ndata, 
                  'redchi': redchi, 'success': success,
                 'ndata' : ndata, 'nvarys' : nvarys,
                 'nfev': nfev}
        
        tabInfo.update(tabFit)

        #write3GaussFitTable(tabInfo,outTable)


    
    def writeOutputTable(self,lineInfo,singleVorBinInfo,result):

        outTableName = self.cfg_par['general']['runNameDir']+'gPlayOut.fits'

        #f os.path.exists(outTableName):
        #    self.updateOutputTable(lineInfo,singleVorBinInfo,tabFit,fitRes,tabErr)
        #else:
        hdr = fits.Header()
        hdr['COMMENT'] = "Here are the outputs of gPlay"
        hdr['COMMENT'] = "Ext 1 = binInfo Ext 2 = fit result Ext 3 = line parameters"
        empty_primary = fits.PrimaryHDU(header=hdr)

        #make extension 1 : vorBinInfo & fitResults
        col1 = fits.Column(name='ID', format='I', array=singleVorBinInfo['ID'])
        col2 = fits.Column(name='BIN_ID', format='I', array=singleVorBinInfo['BIN_ID'])
        col3 = fits.Column(name='X', format='D', array=np.round(singleVorBinInfo['X'],1))
        col4 = fits.Column(name='Y', format='D', array=np.round(singleVorBinInfo['Y'],1))
        col5 = fits.Column(name='X-Pix', format='I', array=singleVorBinInfo['PixX'])
        col6 = fits.Column(name='Y-Pix', format='I', array=singleVorBinInfo['PixX'])


        cols1 = fits.ColDefs([col1, col2, col3, col4, col5, col6])
        
        t1 = fits.BinTableHDU.from_columns(cols1)                
        
        hdul = fits.HDUList([empty_primary,t1])        
        

        aic = result.aic
        bic = result.bic
        redchi = result.redchi
        success = result.success
        ndata = result.ndata
        nvarys = result.nvarys
        nfev = result.nfev
        success = result.success
        
        tabFit = {'aic':aic,'bic': bic, 'ndata': ndata,
          'redchi': redchi, 'success': success,
          'ndata' : ndata, 'nvarys' : nvarys,
          'nfev': nfev}

        modNameList = self.cfg_par['gFit']['modName']
        col7 = fits.Column(name='modName', format='10A', array=modNameList)
        col77 = fits.Column(name='BIN_ID', format='I', array=singleVorBinInfo['BIN_ID'])
        col8 = fits.Column(name='fitSuccess', dtype='L', array=tabFit['success'])
        col9 = fits.Column(name='redChi', format='D', array=tabFit['redchi'])
        col10 = fits.Column(name='aic', format='D', array=tabFit['aic'])
        col11 = fits.Column(name='bic', format='D', array=tabFit['bic'])
        col12 = fits.Column(name='nData', dformat='I', array=tabFit['ndata'])
        col13 = fits.Column(name='nVariables', format='I', array=tabFit['nvarys'])
        col14 = fits.Column(name='nFev', dformat='I', array=tabFit['nfev'])

        cols2 = fits.ColDefs([col7, col77, col8, col9, col10,
            col11, col12, col13])

        t2 = fits.BinTableHDU.from_columns(cols2)
        hdul.append(t2)  

        fitRes = result.params.valuesdict()


        for i in xrange(0,len(lineInfo['ID'])):

            modName = self.cfg_par['gFit']['modName']

            amp_err = result.params[modName+'ln'+str(i)+'_amplitude'].stderr
            sig_err = result.params[modName+'ln'+str(i)+'_sigma'].stderr
            cen_err = result.params[modName+'ln'+str(i)+'_center'].stderr            

            tabErr = {modName+'_amplitude':amp, modName+'_center': ctr, modName+'_sigma': sig,
                  modName+'_fwhm' : fwhm, modName+'_height' : height,
                  modName+'_centreErr':cen_err, modName+'_ampErr': amp_err, 
                  modName+'_sigmaErr': sig_err}


            cL1 = fits.Column(name='Name', format='10A', array=lineInfo['Name'][i])
            cL2 = fits.Column(name='Wave', format='D', array=lineInfo['Wave'][i])
            cL3 = fits.Column(name='g1_Amp', format='D', array=fitRes['g1'+'ln'+str(i)+'_amplitude'])
            cL4 = fits.Column(name='g1_Amp_Err', format='D', array=tabErr['g1_ampErr'])
            cL5 = fits.Column(name='g1_Height', format='D', array=fitRes['g1'+'ln'+str(i)+'_height'])


            g1Ctr = self.lambdaVRad(np.exp(fitRes['g1'+'ln'+str(i)+'_centre']),lineInfo['Wave'][i])
            g1CtrErr = self.lambdaVRad(np.exp(tabErr['g1_centreErr']),lineInfo['Wave'][i])
            g1Sigma = self.lambdaVRad(np.exp(fitRes['g1'+'ln'+str(i)+'_sigma']),lineInfo['Wave'][i])
            g1SigmaErr = self.lambdaVRad(np.exp(tabErr['g1'+'ln'+str(i)+'_sigmaErr']),lineInfo['Wave'][i])
            g1FWHM = self.lambdaVRad(np.exp(fitRes['g1'+'ln'+str(i)+'_fwhm']),lineInfo['Wave'][i])
            
            cL6 = fits.Column(name='g1_Centre', format='D', array=g1Ctr, unit='km/s')
            cL7 = fits.Column(name='g1_Centre_Err', format='D', array=g1CtrErr, unit='km/s')
            cL8 = fits.Column(name='g1_Sigma', format='D', array=g1Sigma, unit='km/s')
            cL9 = fits.Column(name='g1_Sigma_Err', format='D', array=g1Sigma, unit='km/s')
            cL10 = fits.Column(name='g1_Sigma_Err', format='D', array=g1Sigma, unit='km/s')

            colsG1 = fits.ColDefs([cL1, cL2, cL3, cL4,
                    cL5, cL6, cL7,cL8,cL9,cL10])           
            
            if modName == 'g1':

                t3 = fits.BinTableHDU.from_columns(cols2)
                hdul.append(t3)  
            
            else: 
                
                cL11 = fits.Column(name='g2_Amp', format='D', array=fitRes['g2'+'ln'+str(i)+'_amplitude'])
                cL12 = fits.Column(name='g2_Amp_Err', format='D', array=tabErr['g2'+'ln'+str(i)+'_ampErr'])
                cL13 = fits.Column(name='g2_Height', format='D', array=fitRes['g2'+'ln'+str(i)+'_height'])

                g2Ctr = self.lambdaVRad(np.exp(fitRes['g2'+'ln'+str(i)+'_centre']),lineInfo['Wave'][i])
                g2CtrErr = self.lambdaVRad(np.exp(tabErr['g2'+'ln'+str(i)+'_centreErr']),lineInfo['Wave'][i])
                g2Sigma = self.lambdaVRad(np.exp(fitRes['g2'+'ln'+str(i)+'_sigma']),lineInfo['Wave'][i])
                g2SigmaErr = self.lambdaVRad(np.exp(tabErr['g2'+'ln'+str(i)+'_sigmaErr']),lineInfo['Wave'][i])
                g2FWHM = self.lambdaVRad(np.exp(fitRes['g2'+'ln'+str(i)+'_fwhm']),lineInfo['Wave'][i])
            
                cL14 = fits.Column(name='g2_Centre', format='D', array=g2Ctr, unit='km/s')
                cL15 = fits.Column(name='g2_Centre_Err', format='D', array=g2CtrErr, unit='km/s')
                cL16 = fits.Column(name='g2_Sigma', format='D', array=g2Sigma, unit='km/s')
                cL17 = fits.Column(name='g2_Sigma_Err', format='D', array=g2Sigma, unit='km/s')
                cL18 = fits.Column(name='g2_FWHM', format='D', array=g2FWHM, unit='km/s')

                colsG2 = fits.ColDefs([cL11,cL12,cL13,cL14,cL15,
                        cL16,cL17,cL18])

                if self.cfg_par['gFit']['modName'] == 'g2':
                    t3 = fits.BinTableHDU.from_columns([colsG1,colsG2])
                    hdul.append(t3)  
                
                if self.cfg_par['gFit']['modName'] == 'g3':

                    cL19 = fits.Column(name='g3_Amp', format='D', array=fitRes['g2'+'ln'+str(i)+'_amplitude'])
                    cL20 = fits.Column(name='g3_Amp_Err', format='D', array=tabErr['g2_ampErr'])
                    cL21 = fits.Column(name='g3_Height', format='D', array=fitRes['g2'+'ln'+str(i)+'_height'])

                    g3Ctr = self.lambdaVRad(np.exp(fitRes['g3'+'ln'+str(i)+'_centre']),lineInfo['Wave'][i])
                    g3CtrErr = self.lambdaVRad(np.exp(tabErr['g3_centreErr']),lineInfo['Wave'][i])
                    g3Sigma = self.lambdaVRad(np.exp(fitRes['g3'+'ln'+str(i)+'_sigma']),lineInfo['Wave'][i])
                    g3SigmaErr = self.lambdaVRad(np.exp(tabErr['g3_sigmaErr']),lineInfo['Wave'][i])
                    g3FWHM = self.lambdaVRad(np.exp(fitRes['g3'+'ln'+str(i)+'_fwhm']),lineInfo['Wave'][i])
                
                    cL22 = fits.Column(name='g3_Centre', format='D', array=g3Ctr, unit='km/s')
                    cL23 = fits.Column(name='g3_Centre_Err', format='D', array=g3CtrErr, unit='km/s')
                    cL24 = fits.Column(name='g3_Sigma', format='D', array=g3Sigma, unit='km/s')
                    cL25 = fits.Column(name='g3_Sigma_Err', format='D', array=g3Sigma, unit='km/s')
                    cL26 = fits.Column(name='g3_FWHM', format='D', array=g3FWHM, unit='km/s')

                    colsG3 = fits.ColDefs([cL19,cL20,cL21,cL22, cL23, cL24, cL25])
                    
                    t3 = fits.BinTableHDU.from_columns([colsG1,colsG2,colsG3])
                    hdul.append(t3)  
    
        hdul.writeto(outTableName, overwrite=False)

        return 0

    #def writeOutputTable(self,lineInfo,singleVorBinInfo,tabFit,fitRes,tabErr):



#----------------------#
# rc param initialize
#----------------------#
    def loadRcParamsBig(self):
    
        params = {'figure.figsize'      : '10,10',
          'font.family'         :' serif',
          'font.serif'          :'times',
          'font.style'          : 'normal',
          'font.weight'         : 'book',
          'font.size'           : 24,
          'axes.linewidth'      : 2.2,
          'lines.linewidth'     : 2,
          'xtick.labelsize'     : 22,
          'ytick.labelsize'     : 22, 
          'xtick.direction'     :'in',
          'ytick.direction'     :'in',
          'xtick.major.size'    : 6,
          'xtick.major.width'   : 2,
          'xtick.minor.size'    : 3,
          'xtick.minor.width'   : 1,
          'ytick.major.size'    : 6,
          'ytick.major.width'   : 2,
          'ytick.minor.size'    : 3,
          'ytick.minor.width'   : 1, 
          'text.usetex'         : True,
          'text.latex.unicode'  : True
           }
        
        return params

    def loadRcParamsZoom(self):
    
        params = {'figure.figsize'      : '10,10',
          'font.family'         :' serif',
          'font.serif'          :'times',
          'font.style'          : 'normal',
          'font.weight'         : 'book',
          'font.size'           : 13.5,
          'axes.linewidth'      : 2,
          'lines.linewidth'     : 1.8,
          'xtick.labelsize'     : 9,
          'ytick.labelsize'     : 9, 
          'xtick.direction'     :'in',
          'ytick.direction'     :'in',
          'xtick.major.size'    : 4,
          'xtick.major.width'   : 1.5,
          'xtick.minor.size'    : 2,
          'xtick.minor.width'   : 0.75,
          'ytick.major.size'    : 4,
          'ytick.major.width'   : 1.5,
          'ytick.minor.size'    : 2,
          'ytick.minor.width'   : 0.75, 
          'text.usetex'         : True,
          'text.latex.unicode'  : True,
          'legend.fontsize'     : 10
           }
        
        return params

    def plotSpecFit(self,vel,y,result,noise,xx,yy,lineInfo):
        
        velPlot = np.exp(vel)
        yBFit = result.best_fit
        yRes = result.residual
        yInFit = result.init_fit
        key = 'general'
        
        outPlotDir = self.cfg_par[key]['runNameDir']+'/'+self.cfg_par['gPlot']['outPlotDirName']
        if not os.path.exists(outPlotDir):
            os.mkdir(outPlotDir)


        outPlot = outPlotDir+str(xx)+'_'+str(yy)+'_'+self.cfg_par['gFit']['modName']+'.png'       
        

        #dely = result.eval_uncertainty(sigma=3)
            
         # initialize figure
        params = self.loadRcParamsBig()
        plt.rcParams.update(params)
        fig = plt.figure(figsize =(10,8))
        fig.subplots_adjust(hspace=0.0)
        gs = gridspec.GridSpec(1, 1)
        plt.rc('xtick')

        # Initialize subplots
        ax1 = fig.add_subplot(gs[0])

        divider = make_axes_locatable(ax1)
        ax2 = divider.append_axes("bottom", size='15%', pad=0)
        ax1.figure.add_axes(ax2)
        

        #ax.set_xticks([])

        ax1.set_xlabel(r'Wavelength [\AA]')
        ax1.set_ylabel(r'Flux [-]')


        # Calculate axis limits and aspect ratio
        x_min = np.min(velPlot)
        x_max = np.max(velPlot)
        if self.cfg_par['gPlot']['fixed_scale']:
            y1_min = np.min(y)*1.5
            y1_max = np.max(y)*1.5
        else:
            y1_min = np.nanmin(y)*1.3
            y1_max = np.nanmax(y)*1.3

        # Set axis limits
        ax1.set_xlim(x_min, x_max)
        ax1.set_ylim(y1_min, y1_max)

        ax1.tick_params(axis='both', which='major', pad=5)
        #ax1.xaxis.set_minor_locator()
        #ax1.yaxis.set_minor_locator()      
        
        ax1.step(velPlot, y, where='mid', color='black', linestyle='-')
        ax1.plot(velPlot, yBFit, 'r-', label='best fit')
        #ax1.step(vel, yInFit, 'b-', label='init fit')

        #ax1.fill_between(vel, yBFit-dely, yBFit+dely, color="#ABABAB",
        #        label='3-$\sigma$ uncertainty band')
        if self.cfg_par['gFit']['modName'] !='g1':
            comps = result.eval_components()
            for i in xrange(0,len(lineInfo['ID'])):
                
                ax1.plot(velPlot, comps['g1ln'+str(i)+'_'], 'g--')
            
                if self.cfg_par['gFit']['modName'] =='g2':
                    ax1.plot(velPlot, comps['g2ln'+str(i)+'_'], 'm--')    
            
                elif self.cfg_par['gFit']['modName'] !='g2':
                    ax1.plot(velPlot, comps['g2ln'+str(i)+'_'], 'm--')    
                    ax1.plot(velPlot, comps['g3ln'+str(i)+'_'], 'c--')    



        # Calculate axis limits and aspect ratio
        x_min = np.min(velPlot)
        x_max = np.max(velPlot)
        if self.cfg_par['gPlot']['Res-fixed_scale']:
            y1_min = -200
            y1_max = 200
        else:
            y1_min = np.nanmin(yRes)*1.3
            y1_max = np.nanmax(yRes)*1.3

        # Set axis limits
        ax2.set_xlim(x_min, x_max)
        ax2.set_ylim(y1_min, y1_max) 

        #ax2.plot(vel, amp(x, p(x,m,n)))
        ax2.step(velPlot, yRes, 'g-', label='residuals')
        ax2.axhline(color='k', linestyle=':', zorder=0)                           
        ax2.fill_between(velPlot, -noise, noise,
                         facecolor='grey', alpha=0.5,step='mid')

        #ax1.legend(loc='best')



        plt.savefig(outPlot,
                    format='png') # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
        plt.show()
        plt.close()
           
        return 0

    def plotLineZoom(self,vel,y,result,noise,i,j,lineInfo):

        velExp = np.exp(vel)
        yBFit = result.best_fit
        yRes = result.residual
        yInFit = result.init_fit
        key = 'general'
        
        outPlotDir = self.cfg_par[key]['runNameDir']+self.cfg_par['gPlot']['outPlotDirName']
        if not os.path.exists(outPlotDir):
            os.mkdir(outPlotDir)
        
        outPlot = outPlotDir+str(i)+'_'+str(j)+'_SNR-lines.png'       

        params = self.loadRcParamsZoom()
        plt.rcParams.update(params)
        # add one row for the plot with full channel width
        n_plots = len(lineInfo['Wave'])
        n_rows = int(np.ceil(n_plots/3.))
        fig, ax = plt.subplots(squeeze=False,
            ncols=3, nrows=n_rows, figsize=(8.25, 11.67))
        fig.subplots_adjust(hspace=0.)

        #for plot_count in range(n_plots):
        k=0
        for i in xrange(0,len(lineInfo['Wave'])):


            if i == 0:
                j = 0
            elif i % 3 == 0:
                j +=1 
                k =0

            waveMin =  np.log(lineInfo['Wave'][i] - lineInfo['lineRangeAng'][i])
            waveMax =  np.log(lineInfo['Wave'][i] + lineInfo['lineRangeAng'][i])
            idxMin = int(np.where(abs(vel-waveMin)==abs(vel-waveMin).min())[0]) 
            idxMax = int(np.where(abs(vel-waveMax)==abs(vel-waveMax).min())[0] )
        
            x_data_plot = velExp[idxMin:idxMax]
            x_data_plot = self.lambdaVRad(x_data_plot,lineInfo['Wave'][i])
            y_data_plot = y[idxMin:idxMax]
            y_BFit_plot = result.best_fit[idxMin:idxMax]
            y_Res_plot = result.residual[idxMin:idxMax]
            y_sigma_plot = noise[idxMin:idxMax]

            # Calculate axis limits and aspect ratio
            x_min = np.min(x_data_plot)
            x_max = np.max(x_data_plot)
            if self.cfg_par['gPlot']['fixed_scale']:
                y1_min = -75
                y1_max = np.nanmax(y_data_plot)*1.1
            else:
                y1_min = np.nanmin(y_data_plot)*1.1
                y1_max = np.nanmax(y_data_plot)*1.1

            minLoc = np.floor(np.min(x_data_plot)/np.float(self.cfg_par['gPlot']['deltaVlabel']))
            maxLoc = np.ceil(np.max(x_data_plot)/np.float(self.cfg_par['gPlot']['deltaVlabel']))
            xTicks = np.arange(minLoc*self.cfg_par['gPlot']['deltaVlabel'],(maxLoc+1)*self.cfg_par['gPlot']['deltaVlabel'],
                self.cfg_par['gPlot']['deltaVlabel'])
            xTicksStr = [str(hh) for hh in xTicks]

            ax[j][k].set_xticks(xTicks)
            ax[j][k].set_xticklabels([])

            ax[j][k].set_xlim(x_min, x_max)
            ax[j][k].set_ylim(y1_min, y1_max)
            ax[j][k].xaxis.labelpad = 6
            ax[j][k].yaxis.labelpad = 10
            ax[j][k].minorticks_on()
            ax[j][k].tick_params(axis='both', bottom='on', top='on',
                                   left='on', right='on', which='major', direction='in')
            ax[j][k].tick_params(axis='both', bottom='on', top='on',
                                   left='on', right='on', which='minor', direction='in')
            if k==0:
                ylabh = ax[j][k].set_ylabel(
                    r'Flux [-]')
                ylabh.set_verticalalignment('center')

            ax[j][k].tick_params(axis='both', which='major', pad=5)
            #ax1.xaxis.set_minor_locator()
            #ax1.yaxis.set_minor_locator()      
            if lineInfo['Name'][i] == 'Hb':
                lineInfoName =r'H$_\beta$'
            elif lineInfo['Name'][i] == 'Ha':
                lineInfoName =r'H$_\alpha$'
            else:
                lineInfoName =lineInfo['Name'][i]

            ax[j][k].step(x_data_plot, y_data_plot, where='mid', color='black', linestyle='-')
            ax[j][k].plot(x_data_plot, y_BFit_plot, 'r-', label=lineInfoName+str(int(lineInfo['Wave'][i])))
            #ax1.step(vel, yInFit, 'b-', label='init fit')

            #ax1.fill_between(vel, yBFit-dely, yBFit+dely, color="#ABABAB",
            #        label='3-$\sigma$ uncertainty band')
            if self.cfg_par['gFit']['modName'] !='g1':
                comps = result.eval_components()
                for ii in xrange(0,len(lineInfo['ID'])):

                    ax[j][k].plot(x_data_plot, comps['g1ln'+str(ii)+'_'][idxMin:idxMax], 'g--')
                
                    if self.cfg_par['gFit']['modName'] =='g2':
                        ax[j][k].plot(x_data_plot, comps['g2ln'+str(ii)+'_'][idxMin:idxMax], 'm--')    
                
                    elif self.cfg_par['gFit']['modName'] !='g2':
                        ax[j][k].plot(x_data_plot, comps['g2ln'+str(ii)+'_'][idxMin:idxMax], 'm--')    
                        ax[j][k].plot(x_data_plot, comps['g3ln'+str(ii)+'_'][idxMin:idxMax], 'c--')    


            ax[j][k].axvline(color='k', linestyle=':', zorder=0)                           
            ax[j][k].legend(loc='best',handlelength=0.1, handletextpad=0.1,frameon=False)
            ax[j][k].legend(handlelength=0.5)

 
            divider = make_axes_locatable(ax[j][k])
            ax2 = divider.append_axes("bottom", size='20%',pad=0)
            ax2.minorticks_on()

            ax[j][k].figure.add_axes(ax2)
            # Calculate axis limits
            x_min = np.min(x_data_plot)
            x_max = np.max(x_data_plot)
            if self.cfg_par['gPlot']['Res-fixed_scale']:
                y1_min = -200
                y1_max = 200
            else:
                y1_min = np.nanmin(y_Res_plot)*1.3
                y1_max = np.nanmax(y_Res_plot)*1.3

            ax2.set_xticks(xTicks)

            ax2.set_yticks([-150,0,150])
            #ax2.set_yticklabels([])

            # Set axis limits
            ax2.set_xlim(x_min, x_max)
            ax2.set_ylim(y1_min, y1_max) 

            #ax2.plot(vel, amp(x, p(x,m,n)))
            ax2.step(x_data_plot, y_Res_plot, 'g-', label='residuals')
            ax2.axhline(color='k', linestyle=':', zorder=0)                           
            ax2.axvline(color='k', linestyle=':', zorder=0)                           
            ax2.fill_between(x_data_plot, -y_sigma_plot, y_sigma_plot,
                             facecolor='grey', alpha=0.5,step='mid')

           # for the last plot add the x-axis label
            if i == len(lineInfo['Wave'])-1:                
                ax2.set_xticks(xTicks)
                ax2.set_xlabel(
                        r'Velocity\,[$\mathrm{km}\,\mathrm{s}^{-1}$]', labelpad=2)
            else:
                ax2.set_xticklabels([])


            k+=1
        #delete unused subplots
        i+=1
        while not i % 3 ==0:   
            fig.delaxes(ax[j][k])

            ax[j-1][k].axes.get_xaxis().set_ticklabels(xTicksStr)
            ax[j-1][k].axes.get_xaxis().set_ticklabels(xTicksStr)
            k +=1
            i +=1


        plt.savefig(outPlot,dpi=200,bbox_inches='tight',
                    format='png') # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
        plt.show()
        plt.close()                        
