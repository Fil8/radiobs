#!/usr/bin/env python

import sys, os, string
import numpy as np
from astropy.io import fits
import pyregion


import argparse
from  argparse import ArgumentParser
import textwrap as _textwrap

import headPlay
import cvMe

hp = headPlay.headplay()
cv = cvMe.convert()

class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width, initial_indent=indent, subsequent_indent=indent) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text

class fitsplay:

    def maskMe(self,fileName,vals,cutoff,fill_value,region,option):
        
        hh,dd = hp.cleanHead(fileName)


        if option == 'normMode':
            if vals:
                print vals
                xmin=int(vals[0])
                xmax=int(vals[1])
                ymin=int(vals[2])
                ymax=int(vals[3])
                dd[ymin:ymax,xmin:xmax] = 1.0
                index = dd != 1.0
            if cutoff:
                dd[np.isnan(dd)] = cutoff                
                index = dd <= cutoff
                dd[index==False] = 1.0

        if option == 'minorCut':
            if vals:
                xmin=int(vals[0])
                xmax=int(vals[1])
                ymin=int(vals[2])
                ymax=int(vals[3])
                dd[ymin:ymax,xmin:xmax] = 1.0
                index = dd != 1.0
            if cutoff:
                dd[np.isnan(dd)] = cutoff
                index = dd >= cutoff
                dd[index==False] = 1.0

        if option == 'cutInBox':

            xmin=int(vals[0])
            xmax=int(vals[1])
            ymin=int(vals[2])
            ymax=int(vals[3])
            
            dd[0:ymin,:] = fill_value
            dd[ymax:,:] = fill_value
            dd[:,0:xmin] = fill_value
            dd[:,xmax:] = fill_value

            index = dd < cutoff
            dd[index==False] = 1.0

        if option == 'ds9Reg':

            r = pyregion.open(region).as_imagecoord(hh)
            shape = (hh['NAXIS2'], hh['NAXIS1'])
            m = r.get_mask(shape=shape)

            index = m ==False
            dd[m==True] = 1.0

        if option == 'cutInds9':
            
            r = pyregion.open(region).as_imagecoord(hh)
            shape = (hh['NAXIS2'], hh['NAXIS1'])
            m = r.get_mask(shape=shape)

            dd[m==False] = fill_value
            index_cut = dd < cutoff
            
            dd[index_cut==True] = fill_value
            dd[index_cut==False] = 1.0
            index = index_cut

        dd[index==True] = fill_value


        outfile=string.split(fileName,'.fits')[0]
        outfile = outfile+'_mask.fits'

  
        fits.writeto(outfile,dd,hh,overwrite=True)


        return 0

    def centreCut(self,filename,SizeX,SizeY):

        hh,dd = hp.cleanHead(filename)

        xmin = int(np.round(hh['CRPIX1'],0)-np.round(SizeX/2.,0))
        xmax = int(np.round(hh['CRPIX1'],0)+np.round(SizeX/2.,0))
        ymin = int(np.round(hh['CRPIX2'],0)-np.round(SizeY/2.,0))
        ymax = int(np.round(hh['CRPIX2'],0)+np.round(SizeY/2.,0))
        
        newDD = dd[ymin:ymax,xmin:xmax]
        naxis1 = newDD.shape[1]
        naxis2 = newDD.shape[0]
        crpix1 = newDD.shape[1]/2
        crpix2 = newDD.shape[0]/2

        hh['CRPIX1'] = crpix1
        hh['CRPIX2'] = crpix2
        hh['NAXIS1'] = naxis1
        hh['NAXIS2'] = naxis2

        outfile=string.split(filename,'.fits')[0]
        outfile = outfile+'_cutCtr.fits'

        fits.writeto(outfile,newDD,hh,overwrite=True)

        return 0

    def coordCut(self,rap1,decp1,rap2,decp2):

    
        rap1 = cv.hms2deg(rap1)
        rap2 = cv.hms2deg(rap1)
        decp1 = cv.hms2deg(decp1)
        decp2 = cv.hms2deg(decp2)
        hh,dd = hp.cleanHead(fileName)

        w = wcs.WCS(hh)    

        xmin,ymin=w.wcs_world2pix(rap1,decp1,0)
        xmin=int(np.round(xmin,0))
        ymin=int(np.round(ymin,0))
        
        xmax,ymax=w.wcs_world2pix(rap2,decp2,0)
        xmax=int(np.round(xmax,0))
        ymax=int(np.round(ymax,0))
        naxis1=xmax-xmin
        naxis2=ymax-ymin
        
        hh['NAXIS1']=naxis1
        hh['NAXIS2']=naxis2
        hh['CRPIX1']=naxis1/2
        hh['CRPIX2']=naxis2/2
        
        aaa = string.split(filename, '.fits')
        output=aaa[0]+'_coordCut.fits'
        fits.writeto(output,d[ymin:ymax,xmin:xmax],hh,overwrite=True)

        return 0

    def maskCut(self,fileName,maskName):

        hh,dd = hp.cleanHead(fileName)

        mh,mm = hp.cleanHead(maskName)

        index = mm > 0.
        
        dd [index==False] = np.nan

        aaa = string.split(fileName, '.fits')
        output=aaa[0]+'_maskCut.fits'
        fits.writeto(output,dd,hh,overwrite=True)

        return 0

    def multVal(self,fileName,value):

        hh,dd = hp.cleanHead(fileName)

        dd = np.multiply(dd,value)

        aaa = string.split(fileName, '.fits')
        output=aaa[0]+'_mult.fits'
        fits.writeto(output,dd,hh,overwrite=True)

        return 0

    def subFits(self,fileName,subName):

        hh,dd = hp.cleanHead(fileName)

        subh,sub = hp.cleanHead(subName)

        dd = np.subtract(dd,sub)

        aaa = string.split(fileName, '.fits')
        output=aaa[0]+'_sub.fits'
        fits.writeto(output,dd,hh,overwrite=True)

        return 0

    def sumFits(self,fileName,subName):

        hh,dd = hp.cleanHead(fileName)

        subh,sub = hp.cleanHead(subName)

        dd = np.add(dd,sub)

        aaa = string.split(fileName, '.fits')
        output=aaa[0]+'_sum.fits'
        fits.writeto(output,dd,hh,overwrite=True)

        return 0

    def multFits(self,fileName,subName):

        hh,dd = hp.cleanHead(fileName)

        subh,sub = hp.cleanHead(subName)

        dd = np.multiply(dd,sub)

        aaa = string.split(fileName, '.fits')
        output=aaa[0]+'_mult.fits'
        fits.writeto(output,dd,hh,overwrite=True)

        return 0

    def divFits(self,fileName,subName):

        hh,dd = hp.cleanHead(fileName)

        subh,sub = hp.cleanHead(subName)

        dd[dd==0.0] = np.nan

        dd = np.divide(dd,sub)
        dd[np.isinf(dd)] = np.nan

        aaa = string.split(fileName, '.fits')
        output=aaa[0]+'_div.fits'
        fits.writeto(output,dd,hh,overwrite=True)

        return 0

    def sqrtFits(self,fileName):

        hh,dd = hp.cleanHead(fileName)


        dd = np.sqrt(dd)

        aaa = string.split(fileName, '.fits')
        output=aaa[0]+'_sqrt.fits'
        fits.writeto(output,dd,hh,overwrite=True)

        return 0

    def main(self,argv):
        

        for i, arg in enumerate(argv):
            if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

        parser = ArgumentParser(description='radiobs: prtHead print header of fits file',
                                formatter_class=MultilineFormatter,
                                add_help=False)

        add = parser.add_argument

        add("-h", "--help",  action="store_true",
                help="Print help message and exit")

        add("-mm", "--maskMe",  action="store_true",
                help="mask Mask")

        add("-ctrCut", "--centreCut",  action="store_true",
                help="cut subregion from centre of image")

        add("-coCut", "--coordCut",  action="store_true",
                help="cut subregion giving coordinates of low left and up right corner (hms,dms)")

        add("-mCut", "--maskCut",  action="store_true",
                help="cut file based on mask of 1 and 0")

        add("-mVal", "--multVal",  action="store_true",
                help="multiply data by constant")

        add("-subFit", "--subtractFits",  action="store_true",
                help="subtract one image from another")

        add("-sumFit", "--sumFits",  action="store_true",
                help="sum one image by another")

        add("-multFit", "--multiplyFits",  action="store_true",
                help="multiply one image by another")

        add("-divFit", "--divideFits",  action="store_true",
                help="divide one image by another")

        add("-sqrtFit", "--sqrtFits",  action="store_true",
                help="square root of a fits file")

        add('-i', '--input',
            type=str,
            default=False,
            help='''input .fits file''')

        add('-c', '--cutoff',
            type=float,
            default=False,
            help='cutoff value above which select data')

        add('-cInB', '--cutInBox',
            action = 'store_true',
            default=False,
            help='mask region with cutoff within box')

        add('-ds9', '--ds9Reg',
            action = 'store_true',
            default=False,
            help='mask image based on ds9 region')

        add('-cInds9', '--cutInds9',
            action = 'store_true',
            default=False,
            help='mask image with cutoff within ds9 region')

        add('-minCut', '--minorCut',
            action = 'store_true',
            default=False,
            help='mask image values below cutoff')

        add('-val', '--values',
            type=int,
            default=False,
            nargs = 4,
            help='xmin xmax ymin ymax [edges of regions where to cut image in pixels]')

        add('-fVal', '--fillValue',
            type=str,
            default='0',
            help='fill mask with zeros or nan')

        add('-mK', '--multiplyConst',
            type=float,
            default=False,
            help='multiply fits file by this constant')

        add('-x', '--xSize',
            type=int,
            default=False,
            help='xSize of cutout from centre')

        add('-y', '--ySize',
            type=int,
            default=False,
            help='ySize of cutout from centre')

        add('-r', '--region',
            type=str,
            default=False,
            help='ds9 region where to mask image')

        add('-p1', '--lowLeft',
            type=str,
            default=False,
            nargs = 2,
            help='coordinates of low left corner')

        add('-p2', '--upRight',
            type=str,
            default=False,
            nargs = 2,
            help='coordinates of upper right corner')

        add('-mask', '--mask',
            type=str,
            default=False,
            help='''mask .fits file''')

        add('-opFits', '--secondInput',
            type=str,
            default=False,
            help='''second .fits file as input''')

        args, uknown = parser.parse_known_args()

        if args.help:
            print('\n\t************* --- radiobs : fitsPlay : Help --- **************\n')

            parser.print_help()

            print ("""\nRun the following command:\n
        radiobs -fp -mm       -i inputFile.fits -c <cutoff_value> and/or -val <xmin xmax ymin ymax>
        radiobs -fp -mm       -i inputFile.fits  -cutInBox -c <cutoff_value> -val <xmin xmax ymin ymax>
        radiobs -fp -mm       -i inputFile.fits  -ds9 -r <region_name> 
        radiobs -fp -mm       -i inputFile.fits  -cutInds9 -r <region_name> - c <cutoff_value>         
        radiobs -fp -ctrCut   -i inputFile.fits -x <xSize> -y <ySize>
        radiobs -fp -coordCut -i inputFile.fits -p1 <ra dec low left corner (hms)> -p2 <ra dec up right corner (dms)>
        radiobs -fp -maskCut  -i inputFile.fits -mask maskFile.fits
        radiobs -fp -mVal    -i inputFile.fits -mK <constant_value>
        radiobs -fp -subFit  -i inputFile.fits -opFits fileToSubtract.fits
        radiobs -fp -sumFit  -i inputFile.fits -opFits fileToSum.fits
        radiobs -fp -multFit -i inputFile.fits -opFits fileToMultiply.fits
        radiobs -fp -divFit  -i inputFile.fits -opFits fileToDivide.fits
        radiobs -fp -sqrtFit -i inputFile.fits -opFits fileToSquareRoot.fits
                """)
            print('\n\t************* --- radiobs : fitsPlay : DONE --- **************\n')

            sys.exit(0)

        elif args.maskMe:
            print('\n\t************* ---     radiobs : maskMe    --- **************')

            filename = args.input
            vals = args.values
            cutoff = args.cutoff
            region = args.region
            print region
            if args.fillValue == 'nan':
                fillValue = np.nan
            else:
                fillValue = 0

            if args.cutInBox == True:
                option = 'cutInBox'
            elif args.ds9Reg == True:
                option = 'ds9Reg'
            elif args.cutInds9 == True:
                option = 'cutInds9'
            elif args.minorCut == True:
                option = 'minorCut'
            else:
                option = 'normMode'

            self.maskMe(filename,vals,cutoff,fillValue,region,option)
            
            print('\t************* --- radiobs : maskMe : DONE --- **************\n')
        
        elif args.centreCut:
            print('\n\t************* ---     radiobs : ctrCut    --- **************')

            filename = args.input
            x = args.xSize
            y = args.ySize
            self.centreCut(filename,x,y)

            print('\t************* --- radiobs : ctrCut : DONE --- **************\n')

        elif args.coordCut:
            print('\n\t************* ---     radiobs : coordCut    --- **************')
        
            filename = args.input

            rightAscensionP1 = args.p1[0]
            declinationP1 = args.p1[1]
            rightAscensionP1 = rightAscension.replace(',','')

            rightAscensionP2 = args.p2[0]
            declinationP2 = args.p2[1]
            rightAscensionP2 = rightAscension.replace(',','')

            self.centreCut(filename,rightAscensionP1,declinationP1,rightAscensionP2,declinationP2)

            print('\n\t************* --- radiobs : coordCut  : DONE --- **************')

        elif args.maskCut:
            print('\n\t************* ---     radiobs : maskCut    --- **************')
        
            filename = args.input
            maskname = args.mask
            self.maskCut(filename,maskname)

            print('\n\t************* --- radiobs : maskCut  : DONE --- **************')

        elif args.multVal:
            print('\n\t************* ---     radiobs : multVal   --- **************')
        
            filename = args.input
            value = args.multiplyConst
            self.multVal(filename,value)

            print('\n\t************* --- radiobs : multVal  : DONE --- **************')


        elif args.subtractFits:
            print('\n\t************* ---     radiobs : subFits   --- **************')
        
            filename = args.input
            subname = args.secondInput
            self.subFits(filename,subname)

            print('\n\t************* --- radiobs : subFits  : DONE --- **************')

        elif args.multiplyFits:
            print('\n\t************* ---     radiobs : multFits   --- **************')
        
            filename = args.input
            subname = args.secondInput
            self.multFits(filename,subname)

            print('\n\t************* --- radiobs : multFits  : DONE --- **************')


        elif args.sumFits:
            print('\n\t************* ---     radiobs : sumFits   --- **************')
        
            filename = args.input
            subname = args.secondInput
            self.sumFits(filename,subname)

            print('\n\t************* --- radiobs : sumFits  : DONE --- **************')            

        elif args.divideFits:
            print('\n\t************* ---     radiobs : divFits   --- **************')
        
            filename = args.input
            subname = args.secondInput
            self.divFits(filename,subname)

            print('\n\t************* --- radiobs : divFits  : DONE --- **************')     

        elif args.sqrtFits:
            print('\n\t************* ---     radiobs : sqrtFits   --- **************')
        
            filename = args.input
            self.sqrtFits(filename)

            print('\n\t************* --- radiobs : sqrtFits  : DONE --- **************')     


        else:
            print ("""\n\t --- ERROR: Please specify input fits file and command ---""")
            print('\n\t************* --- radiobs : fitsPlay : DONE --- **************\n')
            sys.exit(0)

        print('\n\t************* --- radiobs : fitsPlay : DONE --- **************\n')