#!/usr/bin/env python

import sys, os, string
import numpy as np
from astropy.io import fits

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

class headplay:

    def __init__(self):

        self.rootdir = os.getcwd()+'/'
        self.C = 2.99792458e8 #m/s


    def printHead(self,filename):
        """Function to print the header of a fits file to terminal.
        Parameter:
            filename : str
                Name of input fits file
        Return:
            header : astropu.fits object
                Header of input file
        
        Output:
            Prints header of a fits file to terminal.
        """


        files=fits.open(filename)
        heads=files[0].header

        for i in heads.keys():
            if i != 'HISTORY':
                print(str(i)+'\t'+str(heads[i]))

        files.close()

        return heads

    def printHist(self,filename):

        files=fits.open(filename)

        heads=files[0].header

        if heads['HISTORY']:
            print filename+'\n'
            print heads['HISTORY']
        else:
            print('\n\t --- WARNING: History is not in header ---')

        return 0

    def cleanHead(self,fileName,writeFile):

        base = fits.open(fileName)
        heads = base[0].header
        
        datas = base[0].data
        datas = np.squeeze(datas)

        if 'NAXIS3' in heads:
            del heads['NAXIS3']
        if 'CRVAL3' in heads:
            heads['FREQ'] = heads['CRVAL3']
            del heads['CRVAL3']
        if 'CDELT3' in heads:
            del heads['CDELT3']
        if 'CRPIX3' in heads: 
            del heads['CRPIX3']
        if 'CTYPE3' in heads:        
            del heads['CTYPE3']  
        if 'CROTA3' in heads:
            del heads['CROTA3']
        if 'CUNIT3' in heads:
            del heads['CUNIT3']        
            
        if 'NAXIS4' in heads:
            del heads['NAXIS4']     
        if 'CRVAL4' in heads:        
            del heads['CRVAL4']
        if 'CDELT4' in heads:    
            del heads['CDELT4']
        if 'CRPIX4' in heads:    
            del heads['CRPIX4']
        if 'CTYPE4' in heads:    
            del heads['CTYPE4'] 
        if 'CROTA4' in heads:
            del heads['CROTA4']  
        if 'CUNIT4' in heads:
            del heads['CUNIT4']

        if 'DATAMIN' in heads:
            heads.set('DATAMIN', heads['DATAMIN'], after='BTYPE')
            print 'culo'
        if 'DATAMAC' in heads:    
            heads.set('DATAMAX', heads['DATAMAX'], before='DATAMIN')
        
        if 'CELLSCAL' in heads:        
            heads.set('CELLSCAL', heads['CELLSCAL'], after='DATAMIN')

        if 'RESTFREQ' in heads:
            heads.set('RESTFREQ', heads['RESTFREQ'], after='DATAMIN')
        
        if 'SPECSYS3' in heads:
            heads.set('SPECSYS3', heads['SPECSYS3'], after='RESTFREQ')

        
        if 'WCSAXES' in heads:
            heads['WCSAXES'] = 2
        
        if 'PC1_1' in heads:   
            del heads['PC1_1']            
        if 'PC2_1' in heads:   
            del heads['PC2_1']           
        if 'PC3_1' in heads:   
            del heads['PC3_1']
        if 'PC4_1' in heads:   
            del heads['PC4_1']    
        if 'PC1_2' in heads:
            del heads['PC1_2']
        if 'PC2_2' in heads:
            del heads['PC2_2']
        if 'PC3_2' in heads:
            del heads['PC3_2']
        if 'PC4_2' in heads:
            del heads['PC4_2']
        if 'PC1_3' in heads:
            del heads['PC1_3']
        if 'PC2_3' in heads:
            del heads['PC2_3']
        if 'PC3_3' in heads:    
            del heads['PC3_3']
        if 'PC4_3' in heads:    
            del heads['PC4_3']
        if 'PC1_4' in heads:
            del heads['PC1_4']            
        if 'PC2_4' in heads:
            del heads['PC2_4']
        if 'PC3_4' in heads:
            del heads['PC3_4']
        if 'PC4_4' in heads:
            del heads['PC4_4']
        if 'PV2_2' in heads:
            del heads['PV2_1']
        if 'PV2_2' in heads:
            del heads['PV2_2']
        if 'CD1_1' in heads:
            heads['CDELT1'] = heads['CD1_1']
            del heads['CD1_1']

        if 'CD2_2' in heads:
            heads['CDELT2'] = heads['CD2_2']
            del heads['CD2_2']

        if 'CD1_2' in heads:
            del heads['CD1_2']

        if 'CD2_1' in heads:
            del heads['CD2_1']

        heads['NAXIS'] = 2

        if writeFile == True:
            fits.writeto(fileName,datas,heads,overwrite=True)

        return heads, datas

    def cleanHeadCube(self,fileName,writeFile):
       
        base = fits.open(fileName)
        heads = base[0].header
        
        datas = base[0].data
        datas = np.squeeze(datas)    
            
        if 'NAXIS4' in heads:
            del heads['NAXIS4']     
        if 'CRVAL4' in heads:        
            del heads['CRVAL4']
        if 'CDELT4' in heads:    
            del heads['CDELT4']
        if 'CRPIX4' in heads:    
            del heads['CRPIX4']
        if 'CTYPE4' in heads:    
            del heads['CTYPE4'] 
        if 'CROTA4' in heads:
            del heads['CROTA4']  
        if 'CUNIT4' in heads:
            del heads['CUNIT4']
        
        #if 'WCSAXES' in heads:
        #    heads['WCSAXES'] = 2
        
        if 'PC1_1' in heads:   
            del heads['PC1_1']            
        if 'PC2_1' in heads:   
            del heads['PC2_1']           
        if 'PC3_1' in heads:   
            del heads['PC3_1']
        if 'PC4_1' in heads:   
            del heads['PC4_1']    
        if 'PC1_2' in heads:
            del heads['PC1_2']
        if 'PC2_2' in heads:
            del heads['PC2_2']
        if 'PC3_2' in heads:
            del heads['PC3_2']
        if 'PC4_2' in heads:
            del heads['PC4_2']
        if 'PC1_3' in heads:
            del heads['PC1_3']
        if 'PC2_3' in heads:
            del heads['PC2_3']
        if 'PC3_3' in heads:    
            del heads['PC3_3']
        if 'PC4_3' in heads:    
            del heads['PC4_3']
        if 'PC1_4' in heads:
            del heads['PC1_4']            
        if 'PC2_4' in heads:
            del heads['PC2_4']
        if 'PC3_4' in heads:
            del heads['PC3_4']
        if 'PC4_4' in heads:
            del heads['PC4_4']
        if 'PV2_2' in heads:
            del heads['PV2_1']
        if 'PV2_2' in heads:
            del heads['PV2_2']
        if 'CD1_1' in heads:
            heads['CDELT1'] = heads['CD1_1']
            del heads['CD1_1']

        if 'CD2_2' in heads:
            heads['CDELT2'] = heads['CD2_2']
            del heads['CD2_2']

        if 'CD1_2' in heads:
            del heads['CD1_2']

        if 'CD2_1' in heads:
            del heads['CD2_1']

        heads['NAXIS'] = 3

        if writeFile == True:
            fits.writeto(fileName,datas,heads,overwrite=True)
        print heads
        return heads, datas

    def putHead(self,fileName,key,value,output=False):

        files=fits.open(fileName)

        datas=files[0].data
        heads=files[0].header

        heads[key] = value
        if 'DATAMIN' in heads:
            del heads['DATAMIN']
        if 'DATAMAX' in heads:
            del heads['DATAMAX']

        if output==False:
            output=fileName
        print key,value        

        print heads
        fits.writeto(output,datas,heads,overwrite=True)

        return 0 

    def freqToVrad(self,fileName,output=False):
        """Function to convert datacube from frequency in Hz to VRAD.
        Parameter:
            fileName : str
                Name of input datacube .fits file
        Return:
            heads : 
                Header of converted datacube
        
        Output:
            Saves datacube in VRAD units into working directory.
        """

        heads = self.printHead(fileName)
        f= fits.open(fileName)
        datas=f[0].data
        
        vel = (heads['RESTFRQ']-heads['CRVAL3'])/heads['RESTFRQ']
        
        velStep = (heads['RESTFRQ']-(heads['CRVAL3']+heads['CDELT3']))/heads['RESTFRQ']
        velStep = velStep - vel

        heads['CRVAL3'] = vel*self.C
        heads['CDELT3'] = velStep*self.C
        heads['CTYPE3'] = 'VRAD'

        if output==False:
            aaa = string.split(fileName, '.fits')
            output=aaa[0]+'_vrad.fits'

        fits.writeto(output,datas,heads,overwrite=True)     

        return heads

    def to32Bits(self,fileName):

        ff=fits.open(fileName)
        dd=ff[0].data
        hh=ff[0].header
        dd=dd.astype('float32')

        outfile=string.split(fileName,'.fits')[0]
        outfile = outfile+'_bt32.fits'

        fits.writeto(outfile,dd,hh,overwrite=True)

        return 0

    def main(self,argv):
        

        for i, arg in enumerate(argv):
            if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

        parser = ArgumentParser(description='radiobs: headPlay play with header of fits file',
                                formatter_class=MultilineFormatter,
                                add_help=False)

        add = parser.add_argument

        add("-h", "--help",  action="store_true",
                help="Print help message and exit")

        add("-hp", "--headPlay",  action="store_true",
                default=True,
                help="Call class command")

        add("-prtHd", "--printHead",  action="store_true",
                help="Print header of fits file")

        add("-prtHis", "--printHist",  action="store_true",
                help="Print history of fits file")

        add("-clHd", "--cleanHead",  action="store_true",
                help="Clean header for continuum images and wcs coordinates")

        add("-putHd", "--putHead",  action="store_true",
                help="Put new keyword in header of fits file")

        add("-vrad", "--freqToVrad",  action="store_true",
                help="Convert datacube from frequency to radio velocity")

        add("-to32", "--cvTo32bit",  action="store_true",
                help="Convert fits file to 32 bit")

        add('-i', '--input',
            type=str,
            default=False,
            help='''input .fits file''')

        add('-k', '--key',
            type=str,
            default=False,
            help='name of keyword to insert')

        add('-v', '--val',
            type=str,
            default=False,
            help='value of new keyword')

        args, uknown = parser.parse_known_args()
        if args.help:
            print('\n\t************* --- radiobs : headPlay : Help --- **************\n')

            parser.print_help()

            print ("""\nRun the following command:\n
        radiobs -hp -prtHd  -i inputFile.fits 
        radiobs -hp -prtHis -i inputFile.fits 
        radiobs -hp -clHd   -i inputFile.fits 
        radiobs -hp -to32   -i inputFile.fits 
        radiobs -hp -vrad   -i inputFile.fits
        radiobs -hp -putHd  -i inputFile.fits  -k <key> -val <value> 
                """)
            print('\n\t************* --- radiobs : headPlay : DONE --- **************\n')

            sys.exit(0)

        elif args.printHead:
            print('\n\t************* ---     radiobs : printHead    --- **************')
            print args.input
            filename = args.input
            self.printHead(filename)
        
            print('\n\t************* ---       Header is Over       --- **************')
            print('\n\t************* --- radiobs : printHead : DONE --- **************\n')
        
        elif args.printHist:
            print('\n\t************* ---     radiobs : printHist    --- **************')

            filename = args.input
            self.printHist(filename)

            print('\n\t************* ---      History is Over       --- **************')
            print('\n\t************* --- radiobs : printHist : DONE --- **************\n')

        elif args.cleanHead:
            print('\n\t************* ---     radiobs : cleanHead    --- **************')
            filename = args.input
            self.cleanHead(filename,writeFile=True)
            print('\n\t************* ---       Header is clean      --- **************')
            print('\n\t************* --- radiobs : cleanHead : DONE --- **************\n')

        elif args.putHead:
            print('\n\t************* --- radiobs : putHead  --- **************')
            filename = args.input
            self.putHead(filename,args.key,args.val)
            print('\n\t************* ---   Keword set in header   --- **************')
            print('\n\t************* --- radiobs : putHead : DONE --- **************\n')

        elif args.freqToVrad:
            print('\n\t************* ---     radiobs : freqToVrad    --- **************')

            filename = args.input
            self.freqToVrad(filename)
        
            print('\n\t************* ---       Conversion is Done       --- **************')
            print('\n\t************* --- radiobs : freqToVrad : DONE --- **************\n')

        elif args.cvTo32bit:
            print('\n\t************* ---    radiobs : to32Bits     --- **************')
            filename = args.input
            self.to32Bits(filename)
            print('\n\t************* ---      file converted       --- **************')
            print('\n\t************* --- radiobs : to32Bits : DONE --- **************\n')

        else:
            print ("""\n\t --- ERROR: Please specify input fits file and command ---""")
            print('\n\t************* --- radiobs : headPlay : DONE --- **************\n')
            sys.exit(0)

        print('\n\t************* --- radiobs : headPlay : DONE --- **************\n')



