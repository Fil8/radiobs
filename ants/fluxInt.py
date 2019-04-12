#!/usr/bin/env python

import sys, os, string
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table

import pyregion
from prettytable import PrettyTable

import headPlay

hp = headPlay.headplay()

import argparse
from  argparse import ArgumentParser
import textwrap as _textwrap


class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width, initial_indent=indent, subsequent_indent=indent) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text


class fluxint:
    

    def __init__(self):

        self.rootdir = os.getcwd()+'/'
        #self.filename = sys.argv[1]
    
    def readFreq(self,heads,args):
        '''
        reads frequency from terminal or from header of input fits file
        INPUT:
            heads: header of file
            args: arguments parsed from terminal
        OUTPUT:
            frequency of the observation in MHz
        '''        
        if args.freq:
            freq = float(arg[arg.index('-fr')+1])*1e6
        elif 'CRVAL3' in heads:
            freq =  float(heads['CRVAL3'])
        elif 'FREQ' in heads:
            freq = float(heads['FREQ'])

        else: 
            print'''\t --- ERROR: please specify frequency of observation ---'''
            print '\n\t************* --- radiobs : fluxInt : DONE --- **************\n'
            sys.exit(0)

        return freq

    def maskData(self,slavename,basename):
        '''
        mask datas within ds9 region estimate noise and background outside of it
        INPUT:
            datas: array of data
            heads: header of file
            mask: array of mask
        OUTPUT:
            masked_data
            background
            noise
            masked_number_of_pixels
        '''
        
        base = fits.open(basename)[0]
        slave = fits.open(slavename)[0]

        headername = self.rootdir+'tmp.hdr'
        outname = string.split(slavename,'.fits')[0]+'_regr.fits'

        print '\t ... reproject image to mask ...'
        if not (slave.header['NAXIS1'] == slave.header['NAXIS1'] and base.header['NAXIS2'] == slave.header['NAXIS2']):
            base.header['BMAJ'] =  slave.header['BMAJ']       
            base.header['BMIN'] =  slave.header['BMIN']
            if slave.header['FREQ']: 
                base.header['FREQ'] =  slave.header['FREQ']    
            elif slave.header['CRVAL3']: 
                base.header['FREQ'] =  slave.header['CRVAL3']    

            fits.writeto(basename,base.data,base.header,overwrite=True)
            montage.mGetHdr(basename,headername)
            montage.mProject(slavename,outname,headername)
            os.remove(headername)
            slave_regr = fits.open(outname)[0]
        else:
            slave_regr = fits.open(slavename)[0]
        
        index = base.data > 0

        #mean stats
        background = np.nanmean(slave_regr.data[index==False])
        self.noise = np.nanstd(slave_regr.data[index==False])
        
        pixels = np.count_nonzero(base.data[index==True])
        noise = np.multiply(self.noise,np.sqrt(pixels))
        
        datas = slave_regr.data[index==True]
        

        return datas, background, noise, pixels, slave_regr.header
 

    def maskDatReg(self,dat,heads,region,cutoff):
        '''
        mask datas within ds9 region estimate noise and background outside of it
        INPUT:
            datas: array of data
            heads: header of file
            region: ds9 region
        OUTPUT:
            masked_data
            background
            rms
            masked_number_of_pixels
        '''
        datas = dat.copy()
        print 'aaaaaa'
        # set polygonal mask from ds9 region
        r = pyregion.open(region).as_imagecoord(heads)
        shape = (heads['NAXIS2'], heads['NAXIS1'])
        m = r.get_mask(shape=shape)

        #mean stats
        background = np.nanmean(datas[m==False])
        noise = np.nanstd(datas[m==False])
        
        #if cutoff == 0.0:
        self.cutoff = np.nan
        #elif cutoff < 0:
        #    self.cutoff = noise
        #else:
        #    self.cutoff = cutoff
        #    noise=cutoff            

        datas[m==False] = np.nan
        index_cut = datas < self.cutoff
        
        datas[index_cut] = np.nan
        
        mm = datas.copy()
        mm[:,:] = 1.
        self.pixels = np.count_nonzero(mm[m==True])
        
        #noise = np.multiply(noise,np.sqrt(self.pixels))
        
        return datas, background, noise, self.pixels

    def noiseReg(self,datas,heads,region):
        '''
        mask datas within ds9 region estimate noise and background outside of it
        INPUT:
            datas: array of data
            heads: header of file
            region: ds9 region
        OUTPUT:
            background
            rms
            number_of_pixels_in_region
        '''

        # set polygonal mask from ds9 region
        r = pyregion.open(region).as_imagecoord(heads)
        shape = (heads['NAXIS2'], heads['NAXIS1'])
        m = r.get_mask(shape=shape)

        #mean stats
        background = np.nanmean(datas[m==False])
        noise = np.nanstd(datas[m==False]) 

        mm = datas.copy()

        pixels = np.count_nonzero(mm[m==True])
        return background, noise, pixels

    def noiseMultiReg(self,ldata,lhead,region_dir):
        
        r = [f for f in os.listdir(region_dir) if os.path.isfile(os.path.join(region_dir, f))]

        mean_values = np.zeros(len(r))
        flux_values = np.zeros(len(r))
        pixels = 0.
        for i in xrange(0,len(r)):
            region_name = region_dir+r[i]
            #print region_name
            re = pyregion.open(region_name).as_imagecoord(lhead)
            shape = (lhead['NAXIS2'], lhead['NAXIS1'])
            m_noise = re.get_mask(shape=shape)
            mask_tmp = np.copy(ldata)
            #print np.nansum(mask_tmp[m_noise==True])

            #mask_tmp[m==True] = np.nan
            flux_values[i] = np.nansum(mask_tmp[m_noise==True])
            mean_values[i] = np.nanmean(mask_tmp[m_noise==True])
            pixels += np.count_nonzero(mask_tmp[m_noise==True])
        
        noise = np.nanstd(flux_values)
        back = np.nanmean(flux_values)
        pixels = pixels/len(r)

        return back,noise, pixels


    def measFlux(self,datas,heads,errFlux,option):

        fluxSum=np.nansum(datas)
        if option=='RgCv':
            heads['BMAJ'] = 18.5/3600.
            heads['BMIN'] = 9./3600.

        beamArea = 2*np.pi*float(heads['BMAJ'])*3600./2.35482*float(heads['BMIN'])*3600./2.35482

        pixArea = -float(heads['CDELT2']*3600.)*float(heads['CDELT1']*3600.)

        numPixBeam= beamArea/pixArea


        fluxInt = np.divide(fluxSum,numPixBeam)

        return fluxInt,numPixBeam

    def writeTable(self,heads,fluxInt,noise,numPixBeam,freq,errFlux,outTable='integratedFluxes.tbl'):
        
        #noiseInt = np.divide(noise,numPixBeam)
        noiseInt = float(noise)
        fluxErr = np.sqrt(np.power(fluxInt/100.*errFlux,2)+np.power(noiseInt,2))          
        alldata = np.array([freq*1e-6,np.round(fluxInt,8),np.round(noiseInt,8),np.round(fluxErr,6),np.round(float(heads['BMAJ'])*3600.,3),
            np.round(float(heads['BMIN'])*3600.,3),np.round(float(heads['CDELT2'])*3600.,0),np.round(numPixBeam,3),np.round(self.cutoff*1e3,4),np.round(self.pixels,0)])

        columnames = ['Frequency [MHz]','Integrated Flux [Jy]','Noise [Jy]', 'Error [Jy]', 'BeamMaj [arcsec]','BeamMin [arcsec]',
                'PixSize [arcsec]','Beam/pix','Flux cutoff [mJy/beam]','PixInt']
       
        #self.out_table = self.rootdir+outTable

        if os.path.exists(outTable):
            tt = Table.read(outTable, format='ascii')
            tt.add_row(alldata)
        else: 
            tt = Table(alldata,names=columnames,meta={'name': 'Total flux table'})
            
        ascii.write(tt,outTable, overwrite=True)

        #print table
        alldata = np.array([freq*1e-6,np.round(fluxInt,5),np.round(noiseInt,6),np.round(fluxErr,6),np.round(float(heads['BMAJ'])*3600.,3),
            np.round(float(heads['BMIN'])*3600.,3),np.round(float(heads['CDELT2'])*3600.,0),np.round(numPixBeam,3),np.round(self.cutoff*1e3,4),np.round(self.pixels,0)])

        t = PrettyTable(columnames)
        t.add_row(alldata)
        
        return t, fluxErr


    def main(self,argv):

        for i, arg in enumerate(argv):
            if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

        parser = ArgumentParser(description='radiobs: fluxInt : measure flux of source',
                                formatter_class=MultilineFormatter,
                                add_help=False)

        add = parser.add_argument

        add("-h", "--help",  action="store_true",
                help="Print help message and exit")

        add('-i', '--input',
            type=str,
            default=False,
            help='''input .fits file''')

        add('-t', '--table',
            type=str,
            default=False,
            help='table of filenames')

        add('-m', '--mask',
            type=str,
            default=False,
            help='mask where to measure flux in inputFile')

        add('-r', '--region',
            type=str,
            default=False,
            help='ds9 region where to measure flux in inputFile')

        add('-c', '--cutoff',
            type=float,
            default=1e-3,
            help='cutoff value above which measure flux in inputFile')

        add('-fr', '--freq',
            type=float,
            default=None,
            help='statistical error on flux')

        add('-ferr', '--fluxErr',
            type=float,
            default=5,
            help='statistical error on flux')

        args, unknown = parser.parse_known_args()

        if args.help:
            print '\n\t************* --- radiobs : fluxInt : Help --- **************\n'

            parser.print_help()

            print ("""\nRun the following command:\n
        radiobs -fl -i inputFile.fits -r <ds9RegionName>
        radiobs -fl -i inputFile.fits -m <maskName.fits>
        radiobs -fl -i inputFile.fits -t <table.ascii> -r <ds9RegionName>   
        radiobs -fl -i inputFile.fits -t <table.ascii> -m <maskName>   
                """)
            print '\n\t************* --- radiobs : fluxInt : DONE --- **************\n'

            sys.exit(0)


        if args.region and args.table==False:

            print '\t --- Executing File + Region Combo ---'
            fileName = args.input 
            region  = args.region
            cutoff = args.cutoff
            errFlux = args.fluxErr

            heads,datas=hp.cleanHead(fileName)

            maskedData, background, noise, pixels=self.maskDatReg(datas,heads,region,cutoff)

            freq=self.readFreq(heads,args)
            fluxint, number_pix_beam=self.measFlux(maskedData,heads,noise,errFlux)
            t,flErr = self.writeTable(heads,fluxint,noise,number_pix_beam,freq,errFlux)

            print t

        elif args.mask and args.table==False: 

            print '\t --- Executing File + Mask Combo ---'

            fileName = args.input 
            maskName = args.mask

            maskedData, background, noise, pixels, newHeads = self.maskData(fileName,maskName)
            heads=newHeads

            freq=self.readFreq(heads,args)
            
            fluxint, numPixBeam=self.measFlux(maskedData,heads,noise,errFlux)
            
            t,flErr = self.writeTable(heads,fluxint,noise,numPixBeam,freq,errFlux)

            print t

        if args.mask and args.table: 

            print '\t --- Executing FileList + Mask Combo ---' 
            tableName = args.table
            maskName = args.mask

            for i in xrange(0,len(tableFileNames.columns[0])):

                fileName = tableFileNames.columns[0][i]
                maskedData, background, noise, pixels, newHeads = self.maskData(fileName,maskName)
                heads=newHeads
                
                if 'Frequency' in tableFileNames.dtype.names:
                    freq = tableFileNames.columns[1][i]*1e6

                fluxint, numPixBeam =self.measFlux(maskedData,heads,noise,errFlux)

                t,flErr = self.writeTable(heads,fluxint,noise,numPixBeam,freq,errFlux)
                print t

        if args.region and args.table: 
            print '\t --- Executing FileList + Region Combo ---'

            tableName = args.table
            region = args.region
            cutoff = args.cutoff

            for i in xrange(0,len(tableFileNames.columns[0])):

                fileName = tableFileNames.columns[0][i]
                datas,heads = self.openFile(fileName)
                datas=np.squeeze(datas)
                heads=self.cleanHead(heads)
                maskedData, background, noise, pixels=self.maskDatReg(datas,heads,region,cutoff)
                
                if 'Frequency' in tableFileNames.dtype.names:
                    freq = tableFileNames.columns[1][i]*1e6

                fluxint, numPixBeam =self.measFlux(maskedData,heads,noise,errFlux)

                t,flErr = self.writeTable(heads,fluxint,noise,numPixBeam,freq,errFlux)

                print t


        print '\n\t************* --- radiobs : fluxInt : DONE --- **************\n'