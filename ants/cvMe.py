#!/usr/bin/env python

import sys, os, string
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table

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
        
        decDeg =  dd+mm+ss

        return decDeg

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






