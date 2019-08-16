#!/usr/bin/env python

import sys, os, string
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table

from astropy import units
#from astropy.io import fits
from astropy import wcs

import argparse
from  argparse import ArgumentParser
import textwrap as _textwrap

import headPlay

hp = headPlay.headplay()

class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width, initial_indent=indent, subsequent_indent=indent) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text

class convert:

    def __init__(self):

        self.rootdir = os.getcwd()+'/'

    def hms2deg(self,ra_hms):
        
        ra = string.split(ra_hms, ':')

        hh = float(ra[0])*15
        mm = (float(ra[1])/60)*15
        ss = (float(ra[2])/3600)*15
        
        raDeg = hh+mm+ss

        return raDeg


    # DMS -> degrees
    def dms2deg(self,dec_dms):
        dec = string.split(dec_dms, ':')
        
        dd = abs(float(dec[0]))
        mm = float(dec[1])/60
        ss = float(dec[2])/3600
        
        if float(dec[0])>= 0:
            decDeg = dd+mm+ss
        else:
            decDeg = -(dd+mm+ss)
        
        return decDeg

    def coordToPix(self,fileName,ra,dec,verbose=False):
        '''
        
        Module called by abs_ex
        Converts ra,dec of continuum sources
        into pixel coordinates of the datacube
        
        '''

        #I load the WCS coordinate system:
        #open file

        hh,dd = hp.cleanHead(fileName,writeFile=False)      
        w=wcs.WCS(hh)    
        

        pixels=np.zeros([1,2])
        count_out = 0
        count_flag = 0 
        #for i in xrange(0,len(ra)):
        if ra == 'nan':
            pixels[0, 0]= np.nan
            pixels[0, 1]= np.nan
            count_flag +=1
            if verbose == True:
                print '# Source # '+str([i])+ ' is flagged #'
        else:
            ra_deg = self.hms2deg(ra)
            dec_deg = self.dms2deg(dec)
            px,py=w.wcs_world2pix(ra_deg,dec_deg,0)
            if (0 < round(px,0) < hh['NAXIS1'] and
                    0 < round(py,0) < hh['NAXIS2']): 
                pixels[0, 0]= round(px,0)
                pixels[0, 1]= round(py,0)
            else :
                pixels[0, 0]= np.nan
                pixels[0, 1]= np.nan
                count_out +=1
                if verbose == True:
                    print '# Source #  lies outside the fov of the data cube #'

#        print '# Total number of sources: \t'+str(len(ra))
#        print '# Sources below threshold: \t'+str(count_flag)
#        print '# Sources outside f.o.v.:\t'+str(count_out)
#        print '# Sources to analyze: \t\t'+str(len(ra)-count_flag-count_out)+'\n'

        return pixels
    



    def main(self,argv):
        

        for i, arg in enumerate(argv):
            if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

        parser = ArgumentParser(description='radiobs: prtHead print header of fits file',
                                formatter_class=MultilineFormatter,
                                add_help=False)

        add = parser.add_argument

        add("-h", "--help",  action="store_true",
                help="Print help message and exit")

        add('-ra', '--rightAscension',
            type=str,
            default=False,
            help='Right Ascension to convert')

        add('-dec', '--declination',
            type=str,
            default=False,
            help='Declination to convert')

        add('-rd', '--raDec',
            type=str,
            default=False,
            nargs = 2,
            help='Right ascension and Declination to convert')

        args, uknown = parser.parse_known_args(argv)

        if args.help:
            print('\n\t************* --- radiobs : cvMe : Help --- **************\n')

            parser.print_help()

            print ("""\nRun the following command:\n
        radiobs -cv -ra    <rightAscension>
        radiobs -cv -dec   <Declination>
                """)
            print('\n\t************* --- radiobs : cvMe : DONE --- **************\n')

            sys.exit(0)

        elif args.rightAscension :

            ra = self.hms2deg(args.rightAscension)
            print('\tRa (hms) = {:s}\tRa (deg) = {:f}\n'.format(args.rightAscension,ra))

        elif args.declination:

            dec = self.dms2deg(args.declination)
            print('\tDec (dms) = {:s}\tDec (deg) = {:f}\n'.format(args.declination,dec))

        elif args.raDec:
            rightAscension = args.raDec[0]
            rightAscension = rightAscension.replace(',','')

            declination = args.raDec[1]
            ra = self.hms2deg(rightAscension)
            dec = self.dms2deg(declination)
            print('\tRa  (hms) ={:s}\tRa  (deg) = {:f}\n'.format(rightAscension,ra))            
            print('\tDec (dms) ={:s}\tDec (deg) = {:f}\n'.format(declination,dec))






