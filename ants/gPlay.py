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
   
        runDir = self.cfg_par['general']['workdir']+self.cfg_par['general']['runName']+'/'
        if not os.path.exists(runDir):
            os.mkdir(runDir)
        self.cfg_par['general']['runNameDir'] = runDir

        outTableName = self.cfg_par['general']['runNameDir']+'gPlayOut.fits'

        self.cfg_par['general']['outTableName'] = outTableName

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
        lineInfo = self.openLineList()
        diffusion = 1e-5

        #make outputTable
        #outTable = self.makeOutputTable(lineInfo)
        
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

        #define x-axis array
        lambdaMin = np.log(self.cfg_par['gFit']['lambdaMin'])
        lambdaMax = np.log(self.cfg_par['gFit']['lambdaMax'])


        idxMin = int(np.where(abs(wave-lambdaMin)==abs(wave-lambdaMin).min())[0]) 
        idxMax = int(np.where(abs(wave-lambdaMax)==abs(wave-lambdaMax).min())[0] )


        Ydim = dd.shape[1]
        Xdim = dd.shape[2]
        
        binID, binArr, fitResArr, lineArr = self.makeInputArrays(lineInfo, Xdim, Ydim)
       
        counter = 0
        #for j in xrange(205,208):
        #    for i in xrange(250,252):
        #print dd.shape[1]*dd.shape[2]
        #sys.exit(0)
        for j in xrange(0,dd.shape[1]):
            for i in xrange(0,dd.shape[2]):
                #print 'ciao'
                y = dd[idxMin:idxMax,j,i]

                waveCut = wave[idxMin:idxMax]
                #check if spectrum is not empty                   
                if np.sum(y)>0:

                    gMod,gPars = self.lineModDef(waveCut,y,lineInfo)

                    # identify voronoi bin
                    xVal = xAxis[i]
                    yVal = yAxis[j]
                    
                    index = np.where((vorBinInfo['X'] < (xVal+pxSize/2.+diffusion)) & 
                    ((xVal-pxSize/2.-diffusion) < vorBinInfo['X']) & (vorBinInfo['Y'] < (yVal+pxSize/2.+diffusion)) & 
                    ((yVal-pxSize/2.-diffusion) < vorBinInfo['Y']))
                    #print index
                    if np.sum(index)>0: 
                        binArr = self.updateBinArray(binArr,vorBinInfo,index,i,j,counter)
                        binIDName = binArr['BIN_ID'][counter]     
                    else:
                        #print 'continue'
                        fitResArr = np.delete(fitResArr,counter,0)
                        lineArr = np.delete(lineArr,counter,0)  
                        counter+=1
                        continue
                    
                    #check if it is first time in bin
                    if binIDName not in binID[:,:] and np.sum(index)>0:

                        
                        binID[j,i] = binIDName
                        noiseVec = noiseBin[binIDName][:]

                        # FIT
                        result = gMod.fit(y, gPars, x=waveCut)

                        fitResArr = self.updateFitArray(fitResArr,result,binIDName,counter)
                        lineArr = self.updateLineArray(lineArr,result,lineInfo,binIDName,counter)

                        #plot Fit
                        #self.plotSpecFit(waveCut, y,result,noiseVec[idxMin:idxMax],i,j,lineInfo,vorBinInfo[index])
                        #self.plotLineZoom(waveCut, y,result,noiseVec[idxMin:idxMax],i,j,lineInfo,vorBinInfo[index])
                    else:
                        p#rint counter, i,j
                        fitResArr = np.delete(fitResArr,counter,0)
                        lineArr = np.delete(lineArr,counter,0)                                
                else:

                    fitResArr = np.delete(fitResArr,counter,0)
                    lineArr = np.delete(lineArr,counter,0)                                

                counter+=1
                #print 'end_for'
        self.saveOutputTable(binArr, fitResArr, lineArr)
    
        print('''\t+---------+\n\t gFit done\n\t+---------+''')
    
        return 0
    
    def openLineList(self):
        workDir = self.cfg_par['general']['workdir']
       
        lineList = workDir+self.cfg_par['general']['lineListName']
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

    def makeInputArrays(self,lineInfo, Xdim,Ydim):

        binID = np.zeros([Ydim,Xdim],dtype=int)

        nam = tuple (['ID', 'BIN_ID', 'X', 'Y', 'PixX', 'PixY'])
        binArr = np.empty([Ydim*Xdim], dtype={'names':nam,
                          'formats':('i4', 'i4', 'i4', 'f8', 'f8', 'i4', 'i4')})
        nam = tuple(['BIN_ID', 'fitSuccess', 'redChi', 'aic', 'bic', 'nData', 'nVariables', 'nFev'])
        fitResArr = np.zeros([Ydim*Xdim], dtype={'names':nam,
                          'formats':( 'i4', '?', 'f8', 'f8', 'f8', 'i4', 'i4', 'i4')})

        lineNameList = []
        frmList = []
        for i in xrange (0,len(lineInfo['ID'])):
            lineName = str(lineInfo['Name'][i])+str(int(lineInfo['Wave'][i]))
            if '[' in lineName:
                lineName = lineName.replace("[", "")
                lineName = lineName.replace("]", "")

            lineNameList.append(lineName)
            frmList.append('i4')
            lineNameList.append('g1_Amp_'+lineName)
            frmList.append('f8')
            lineNameList.append('g1_Height_'+lineName)
            frmList.append('f8')
            lineNameList.append('g1_Centre_'+lineName)
            frmList.append('f8')
            lineNameList.append('g1_Sigma_'+lineName)
            frmList.append('f8')
            lineNameList.append('g1_FWHM_'+lineName)
            frmList.append('f8')

            
            if self.cfg_par['gFit']['modName'] == 'g2':
                
                lineNameList.append('g2_Amp_'+lineName)
                frmList.append('f8')
                lineNameList.append('g2_Height_'+lineName)
                frmList.append('f8')
                lineNameList.append('g2_Centre_'+lineName)
                frmList.append('f8')
                lineNameList.append('g2_Sigma_'+lineName)
                frmList.append('f8')
                lineNameList.append('g2_FWHM_'+lineName)
                frmList.append('f8')

            
                if self.cfg_par['gFit']['modName'] == 'g3':

                    lineNameList.append('g3_Amp_'+lineName)
                    frmList.append('f8')
                    lineNameList.append('g3_Height_'+lineName)
                    frmList.append('f8')
                    lineNameList.append('g3_Centre_'+lineName)
                    frmList.append('f8')
                    lineNameList.append('g3_Sigma_'+lineName)
                    frmList.append('f8')
                    lineNameList.append('g3_FWHM_'+lineName)
                    frmList.append('f8')
        
        if self.cfg_par['gFit']['modName'] == 'g1':
            lineArr = np.zeros([Ydim*Xdim], dtype={'names':(lineNameList), 'formats':(frmList)})
        elif self.cfg_par['gFit']['modName'] == 'g2':
            lineArr = np.zeros([Ydim*Xdim], dtype={'names':(lineNameList), 'formats':(frmList)})
        elif self.cfg_par['gFit']['modName'] == 'g3':
            lineArr = np.zeros([Ydim*Xdim], dtype={'names':(lineNameList), 'formats':(frmList)})

        return binID, binArr, fitResArr, lineArr


    def updateBinArray(self,binArr,vorBinInfo,index,i,j,counter):
  
        binArr['BIN_ID'][counter] = vorBinInfo['BIN_ID'][index]
        binArr['ID'][counter] = vorBinInfo['ID'][index]
        binArr['X'][counter] = vorBinInfo['X'][index]
        binArr['Y'][counter] = vorBinInfo['Y'][index]
        binArr['PixX'][counter] = int(i)
        binArr['PixY'][counter] = int(j)

        return binArr


    def updateFitArray(self,fitResArr,result,binIDName,counter):

        aic = result.aic
        bic = result.bic
        redchi = result.redchi
        success = result.success
        ndata = result.ndata
        nvarys = result.nvarys
        nfev = result.nfev
        success = result.success
  
        fitResArr['BIN_ID'][counter] = binIDName
        fitResArr['fitSuccess'][counter] = success
        fitResArr['redChi'][counter] = redchi
        fitResArr['aic'][counter] = aic
        fitResArr['nData'][counter] = bic
        fitResArr['nVariables'][counter] = nvarys
        fitResArr['nFev'][counter] = nfev
        fitResArr['nData'][counter] = ndata
        
        return fitResArr

    def updateLineArray(self,lineArr,result,lineInfo,binIDName,counter):
        
        fitRes = result.params.valuesdict()


        modName = self.cfg_par['gFit']['modName']

        for ii in xrange(0,len(lineInfo['ID'])):

            lineName = str(lineInfo['Name'][ii])
            if '[' in lineName:
                lineName = lineName.replace("[", "")
                lineName = lineName.replace("]", "")
            
            lineName = lineName+str(int(lineInfo['Wave'][ii]))

            amp = fitRes['g1ln'+str(ii)+'_amplitude']
            ctr = fitRes['g1ln'+str(ii)+'_center']
            sig = fitRes['g1ln'+str(ii)+'_sigma']
            fwhm = fitRes['g1ln'+str(ii)+'_fwhm']
            height = fitRes['g1ln'+str(ii)+'_height']

            g1Ctr = self.lambdaVRad(np.exp(ctr),lineInfo['Wave'][ii])
            g1Sigma = self.lambdaVRad(np.exp(sig),lineInfo['Wave'][ii])
            g1FWHM = self.lambdaVRad(np.exp(fwhm),lineInfo['Wave'][ii])

            #amp_err = result.params[modName+'ln'+str(i)+'_amplitude'].stderr
            #sig_err = result.params[modName+'ln'+str(i)+'_sigma'].stderr
            #g1SigmaErr = self.lambdaVRad(np.exp(sig_err),lineInfo['Wave'][i])
            #cen_err = result.params[modName+'ln'+str(i)+'_center'].stderr  
            #g1CtrErr = self.lambdaVRad(np.exp(cen_err),lineInfo['Wave'][i])
            lineArr[lineName][counter] = int(lineInfo['Wave'][ii])     
            lineArr['g1_Amp_'+lineName][counter] = amp
            lineArr['g1_Height_'+lineName][counter] = height
            lineArr['g1_Centre_'+lineName][counter] = g1Ctr
            lineArr['g1_Sigma_'+lineName][counter] = g1Sigma
            lineArr['g1_FWHM_'+lineName][counter] = g1FWHM

            if modName != 'g1':

                amp = fitRes['g2ln'+str(i)+'_amplitude']
                ctr = fitRes['g2ln'+str(i)+'_center']
                sig = fitRes['g2ln'+str(i)+'_sigma']
                fwhm = fitRes['g2ln'+str(i)+'_fwhm']
                height = fitRes['g2ln'+str(i)+'_height']

                g2Ctr = self.lambdaVRad(np.exp(ctr),lineInfo['Wave'][ii])
                g2Sigma = self.lambdaVRad(np.exp(sig),lineInfo['Wave'][ii])
                g2FWHM = self.lambdaVRad(np.exp(fwhm),lineInfo['Wave'][ii])

                #amp_err = result.params[modName+'ln'+str(i)+'_amplitude'].stderr
                #sig_err = result.params[modName+'ln'+str(i)+'_sigma'].stderr
                #g1SigmaErr = self.lambdaVRad(np.exp(sig_err),lineInfo['Wave'][i])
                #cen_err = result.params[modName+'ln'+str(i)+'_center'].stderr  
                #g1CtrErr = self.lambdaVRad(np.exp(cen_err),lineInfo['Wave'][i])
          
                lineArr['g2_Amp_'+lineName][counter] = amp
                lineArr['g2_Height_'+lineName][counter] = height
                lineArr['g2_Centre_'+lineName][counter] = g1Ctr
                lineArr['g2_Sigma_'+lineName][counter] = g1Sigma
                lineArr['g2_FWHM_'+lineName][counter] = g1FWHM

                if modName == 'g3':

                    amp = fitRes['g3ln'+str(i)+'_amplitude']
                    ctr = fitRes['g3ln'+str(i)+'_center']
                    sig = fitRes['g3ln'+str(i)+'_sigma']
                    fwhm = fitRes['g3ln'+str(i)+'_fwhm']
                    height = fitRes['g3ln'+str(i)+'_height']

                    g3Ctr = self.lambdaVRad(np.exp(ctr),lineInfo['Wave'][ii])
                    g3Sigma = self.lambdaVRad(np.exp(sig),lineInfo['Wave'][ii])
                    g3FWHM = self.lambdaVRad(np.exp(fwhm),lineInfo['Wave'][ii])

                    #amp_err = result.params[modName+'ln'+str(i)+'_amplitude'].stderr
                    #sig_err = result.params[modName+'ln'+str(i)+'_sigma'].stderr
                    #g1SigmaErr = self.lambdaVRad(np.exp(sig_err),lineInfo['Wave'][i])
                    #cen_err = result.params[modName+'ln'+str(i)+'_center'].stderr  
                    #g1CtrErr = self.lambdaVRad(np.exp(cen_err),lineInfo['Wave'][i])
              
                    lineArr['g3_Amp_'+lineName][counter] = amp
                    lineArr['g3_Height_'+lineName][counter] = height
                    lineArr['g3_Centre_'+lineName][counter] = g3Ctr
                    lineArr['g3_Sigma_'+lineName][counter] = g3Sigma
                    lineArr['g3_FWHM_'+lineName][counter] = g3FWHM
            

        return lineArr

    def saveOutputTable(self, binArr, fitResArr, lineArr):
        
        outTableName = self.cfg_par['general']['runNameDir']+'/gPlayOut1.fits'
        modNameList = self.cfg_par['gFit']['modName']

        if os.path.exists(self.cfg_par['general']['outTableName']):
            hdul = fits.open(self.cfg_par['general']['outTableName'])
            t2 = fits.BinTableHDU.from_columns(fitResArr,name='FitRes_'+modNameList)
            hdul.append(t2)  
            t3 = fits.BinTableHDU.from_columns(lineArr,name='LineRes_'+modNameList)
            hdul.append(t3)  
        else:    
            hdr = fits.Header()
            hdr['COMMENT'] = "Here are the outputs of gPlay"
            hdr['COMMENT'] = "Ext 1 = binInfo Ext 2 = fit result Ext 3 = line parameters"
            
            empty_primary = fits.PrimaryHDU(header=hdr)
           
            t1 = fits.BinTableHDU.from_columns(binArr,name='BinInfo')  
            hdul = fits.HDUList([empty_primary,t1])        

            t2 = fits.BinTableHDU.from_columns(fitResArr,name='FitRes_'+modNameList)
            hdul.append(t2)  

            t3 = fits.BinTableHDU.from_columns(lineArr,name='LineRes_'+modNameList)
            hdul.append(t3)  

        hdul.writeto(outTableName,overwrite=True)

        return

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
    
        params = {
          'font.family'         :' serif',
          'font.serif'          :'times',
          'font.style'          : 'normal',
          'font.weight'         : 'book',
          'font.size'           : 12,
          'axes.linewidth'      : 2,
          'lines.linewidth'     : 1.5,
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

    def plotSpecFit(self,vel,y,result,noise,xx,yy,lineInfo,singleVorBinInfo):
        
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

        ax1.set_xlabel(r'Wavelength [$\AA$]')
        ax1.set_ylabel(r'Flux [-]')


        # Calculate axis limits and aspect ratio
        x_min = np.min(velPlot)
        x_max = np.max(velPlot)
        if self.cfg_par['gPlot']['fixed_scale']:
            y1_min = np.min(y)*1.2
            y1_max = np.max(y)*1.2
        else:
            y1_min = np.nanmin(y)*1.1
            y1_max = np.nanmax(y)*1.1

        # Set axis limits
        ax1.set_xlim(x_min, x_max)
        ax1.set_ylim(y1_min, y1_max)

        ax1.tick_params(axis='both', which='major', pad=5)
        #ax1.xaxis.set_minor_locator()
        #ax1.yaxis.set_minor_locator()      
        
        ax1.step(velPlot, y, where='mid', color='black', linestyle='-')
        ax1.plot(velPlot, yBFit, 'r-', label='best fit')
        #ax1.step(vel, yInFit, 'b-', label='init fit')
        aicStr = str(int(result.aic))
        bicStr = str(int(result.bic))
        redchiStr = str(int(result.redchi))
        successStr = str(result.success)       
        xText = self.cfg_par['gFit']['lambdaMin']+50


        ax1.text(xText, y1_max*0.90, r'BIN ID:\t'+str(singleVorBinInfo['BIN_ID'][0]), {'color': 'k', 'fontsize': 20})
        ax1.text(xText, y1_max*0.94, r'X,Y:\t'+str(xx)+','+str(yy), {'color': 'k', 'fontsize': 20})

        #ax1.text(xText, x_max*0.85, r'Success:\t'+successStr, {'color': 'b'})
        ax1.text(xText, y1_max*0.88, r'$\tilde{\chi}^2$:\t'+redchiStr, {'color': 'k', 'fontsize': 20})
        ax1.text(xText, y1_max*0.82, r'aic:\t'+aicStr, {'color': 'k', 'fontsize': 20})
        ax1.text(xText, y1_max*0.76, r'bic:\t'+bicStr, {'color': 'k', 'fontsize': 20})

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
            y1_min = np.nanmin([-200.,np.nanmin(-noise)*1.1,np.nanmin(-yRes)*1.1])
            y1_max = np.nanmax([+200.,np.nanmax(+noise)*1.1,np.nanmax(+yRes)*1.1])
        else:
            y1_min = np.nanmin(yRes)*1.1
            y1_max = np.nanmax(yRes)*1.1  

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

    def plotLineZoom(self,vel,y,result,noise,i,j,lineInfo,singleVorBinInfo):

        velExp = np.exp(vel)
        yBFit = result.best_fit
        yRes = result.residual
        yInFit = result.init_fit
        key = 'general'
        
        outPlotDir = self.cfg_par[key]['runNameDir']+'/'+self.cfg_par['gPlot']['outPlotDirName']
        if not os.path.exists(outPlotDir):
            os.mkdir(outPlotDir)
        
        outPlot = outPlotDir+str(singleVorBinInfo['BIN_ID'][0])+'_'+self.cfg_par['gFit']['modName']+'.png'       

        params = self.loadRcParamsZoom()
        plt.rcParams.update(params)
        # add one row for the plot with full channel width
        n_plots = len(lineInfo['Wave'])
        n_rows = int(np.ceil(n_plots/3.))
        #fig, ax = plt.subplots(squeeze=False,
        #    ncols=3, nrows=n_rows, figsize=(8.25, 11.67))


        fig = plt.figure(figsize=(8.25, 11.67), constrained_layout=True)
        fig.subplots_adjust(hspace=0.)

        gs_top = plt.GridSpec(nrows=n_rows+1, ncols=3,  figure=fig, top=0.95)
        gs_base = plt.GridSpec(nrows=n_rows+1, ncols=3,  figure=fig, hspace=0,top=0.91)


        #gs = fig.add_gridspec(nrows=n_rows+1, ncols=3, left=0.05, figsize=(8.25, 11.67))
        #gs = gridspec.GridSpec(nrows=n_rows+1, ncols=3,  figure=fig)
        
        wave_ax = self.addFullSubplot(fig,gs_top,vel,y,result,noise,i,j,lineInfo,singleVorBinInfo)
        
 
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
            
            if np.abs(xTicks[0]) > xTicks[-1]:
                    xTicks=np.delete(xTicks,0)
            elif np.abs(xTicks[0]) < xTicks[-1]:
                    xTicks = np.delete(xTicks,-1)


            xTicksStr = [str(hh) for hh in xTicks]
        
            ax = fig.add_subplot(gs_base[j+1,k])
            

            ax.set_xticks(xTicks)
            ax.set_xticklabels([])

            ax.set_xlim(x_min, x_max)
            ax.set_ylim(y1_min, y1_max)
            ax.xaxis.labelpad = 6
            ax.yaxis.labelpad = 10
            ax.minorticks_on()
            ax.tick_params(axis='both', bottom='on', top='on',
                                   left='on', right='on', which='major', direction='in')
            ax.tick_params(axis='both', bottom='on', top='on',
                                   left='on', right='on', which='minor', direction='in')
            if k==0:
                ylabh = ax.set_ylabel(
                    r'Flux [-]')
                ylabh.set_verticalalignment('center')

            ax.tick_params(axis='both', which='major', pad=5)
            #ax1.xaxis.set_minor_locator()
            #ax1.yaxis.set_minor_locator()      
            if lineInfo['Name'][i] == 'Hb':
                lineInfoName =r'H$_\beta$'
            elif lineInfo['Name'][i] == 'Ha':
                lineInfoName =r'H$_\alpha$'
            else:
                lineInfoName =lineInfo['Name'][i]

            ax.step(x_data_plot, y_data_plot, where='mid', color='black', linestyle='-')
            ax.plot(x_data_plot, y_BFit_plot, 'r-', label=lineInfoName+str(int(lineInfo['Wave'][i])))
            #ax2.fill_between(0, y_BFit_plot, y_sigma_plot,
            #                 facecolor='grey', alpha=0.5,step='mid')
            #ax1.step(vel, yInFit, 'b-', label='init fit')

            #ax1.fill_between(vel, yBFit-dely, yBFit+dely, color="#ABABAB",
            #        label='3-$\sigma$ uncertainty band')
            if self.cfg_par['gFit']['modName'] !='g1':
                comps = result.eval_components()
                for ii in xrange(0,len(lineInfo['ID'])):

                    ax.plot(x_data_plot, comps['g1ln'+str(ii)+'_'][idxMin:idxMax], 'g--')
                
                    if self.cfg_par['gFit']['modName'] =='g2':
                        ax.plot(x_data_plot, comps['g2ln'+str(ii)+'_'][idxMin:idxMax], 'm--')    
                
                    elif self.cfg_par['gFit']['modName'] !='g2':
                        ax.plot(x_data_plot, comps['g2ln'+str(ii)+'_'][idxMin:idxMax], 'm--')    
                        ax.plot(x_data_plot, comps['g3ln'+str(ii)+'_'][idxMin:idxMax], 'c--')    


            ax.axvline(color='k', linestyle=':', zorder=0)                           
            legend = ax.legend(loc='best',handlelength=0.0, handletextpad=0.0,frameon=False)
            legend.get_frame().set_facecolor('none')

 
            divider = make_axes_locatable(ax)
            ax2 = divider.append_axes("bottom", size='20%',pad=0)
            ax2.minorticks_on()

            ax.figure.add_axes(ax2)
            # Calculate axis limits
            x_min = np.nanmin(x_data_plot)
            x_max = np.nanmax(x_data_plot)
            if self.cfg_par['gPlot']['Res-fixed_scale']:
                y1_min = np.nanmin([-200.,np.nanmin(-y_sigma_plot)*1.5,np.nanmin(-y_Res_plot)*1.5])
                y1_max = np.nanmax([+200.,np.nanmax(+y_sigma_plot)*1.5,np.nanmax(+y_Res_plot)*1.5])
            else:
                y1_min = np.nanmin(y_Res_plot)*1.1
                y1_max = np.nanmax(y_Res_plot)*1.1

            ax2.set_xticks(xTicks)

            #ax2.set_yticks([-150,0,150])
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
            if i >= len(lineInfo['Wave'])-3:                
                ax2.set_xticks(xTicks)
                ax2.set_xlabel(
                        r'v$_\mathrm{radio}$\,[$\mathrm{km}\,\mathrm{s}^{-1}$]', labelpad=2)
            else:
                ax2.set_xticklabels([])


            k+=1
        #delete unused subplots
        #i+=1
        #while not i % 3 ==0:   
        #    fig.delaxes(ax)

        #    gs[j-1][k].axes.get_xaxis().set_ticklabels(xTicksStr)
        #    gs[j-1][k].axes.get_xaxis().set_ticklabels(xTicksStr)
        #    k +=1
        #    i +=1


        plt.savefig(outPlot,dpi=200,bbox_inches='tight',
                    format='png') # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
        #plt.show()
        plt.close()  

    def addFullSubplot(self,fig,gs,vel,y,result,noise,xx,yy,lineInfo,singleVorBinInfo):
        
        velPlot = np.exp(vel)
        yBFit = result.best_fit
        yRes = result.residual
        yInFit = result.init_fit
        key = 'general'
        
        #outPlotDir = self.cfg_par[key]['runNameDir']+'/'+self.cfg_par['gPlot']['outPlotDirName']
        #if not os.path.exists(outPlotDir):
        #    os.mkdir(outPlotDir)


        #outPlot = outPlotDir+str(xx)+'_'+str(yy)+'_'+self.cfg_par['gFit']['modName']+'.png'       
        

        #dely = result.eval_uncertainty(sigma=3)
            
         # initialize figure
        #params = self.loadRcParamsBig()
        #plt.rcParams.update(params)
        #fig = plt.figure(figsize =(10,8))
        #fig.subplots_adjust(hspace=0.0)
        #gs = gridspec.GridSpec(1, 1)
        #plt.rc('xtick')

        # Initialize subplots
        ax1 = fig.add_subplot(gs[0,:])

        divider = make_axes_locatable(ax1)
        ax2 = divider.append_axes("bottom", size='15%', pad=0)
        ax1.figure.add_axes(ax2)
        

        #ax.set_xticks([])

        ax1.set_xlabel(r'Wavelength [$\AA$]',labelpad=15)
        ax1.set_ylabel(r'Flux [-]')


        # Calculate axis limits and aspect ratio
        x_min = np.min(velPlot)
        x_max = np.max(velPlot)
        if self.cfg_par['gPlot']['fixed_scale']:
            y1_min = np.min(y)*1.2
            y1_max = np.max(y)*1.2
        else:
            y1_min = np.nanmin(y)*1.1
            y1_max = np.nanmax(y)*1.1

        # Set axis limits
        ax1.set_xlim(x_min, x_max)
        ax1.set_ylim(y1_min, y1_max)

        ax1.tick_params(axis='both', which='major', pad=5)
        #ax1.xaxis.set_minor_locator()
        #ax1.yaxis.set_minor_locator()      
        
        ax1.step(velPlot, y, where='mid', color='black', linestyle='-')
        ax1.plot(velPlot, yBFit, 'r-', label='best fit')
        #ax1.step(vel, yInFit, 'b-', label='init fit')
        aicStr = str(int(result.aic))
        bicStr = str(int(result.bic))
        redchiStr = str(int(result.redchi))
        successStr = str(result.success)       
        xText = self.cfg_par['gFit']['lambdaMin']+50


        ax1.text(xText, y1_max*0.90, r'BIN ID:\t'+str(singleVorBinInfo['BIN_ID'][0]), {'color': 'k', 'fontsize': 8})
        ax1.text(xText, y1_max*0.80, r'X,Y:\t'+str(xx)+','+str(yy), {'color': 'k', 'fontsize': 8})

        #ax1.text(xText, x_max*0.85, r'Success:\t'+successStr, {'color': 'b'})
        ax1.text(xText, y1_max*0.70, r'$\tilde{\chi}^2$:\t'+redchiStr, {'color': 'k', 'fontsize': 8})
        ax1.text(xText, y1_max*0.60, r'aic:\t'+aicStr, {'color': 'k', 'fontsize': 8})
        ax1.text(xText, y1_max*0.50, r'bic:\t'+bicStr, {'color': 'k', 'fontsize': 8})

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
            y1_min = np.nanmin([-200.,np.nanmin(-noise)*1.1,np.nanmin(-yRes)*1.1])
            y1_max = np.nanmax([+200.,np.nanmax(+noise)*1.1,np.nanmax(+yRes)*1.1])
        else:
            y1_min = np.nanmin(yRes)*1.1
            y1_max = np.nanmax(yRes)*1.1  

        # Set axis limits
        ax2.set_xlim(x_min, x_max)
        ax2.set_ylim(y1_min, y1_max) 

        #ax2.plot(vel, amp(x, p(x,m,n)))
        ax2.step(velPlot, yRes, 'g-', label='residuals')
        ax2.axhline(color='k', linestyle=':', zorder=0)                           
        ax2.fill_between(velPlot, -noise, noise,
                         facecolor='grey', alpha=0.5,step='mid')

        #ax1.legend(loc='best')



        #plt.savefig(outPlot,
        #            format='png') # if pdf,dpi=300,transparent=True,bbox_inches='tight',overwrite=True)
        #plt.show()
        #plt.close()
           
        return ax1

    def main(self,argv):
        
        for i, arg in enumerate(argv):
            if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

        parser = ArgumentParser(description='radiobs: gPlay play with gaussian fits in datacubes',
                                formatter_class=MultilineFormatter,
                                add_help=False)

        add = parser.add_argument

        add("-h", "--help",  action="store_true",
                help="Print help message and exit")

        print('\n\t************* --- radiobs : gPlay : DONE --- **************\n')

