#!/usr/bin/env python

import sys, os
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

    def printHead(self,filename):

        files=fits.open(filename)
        heads=files[0].header

        for i in heads.keys():
            if i != 'HISTORY':
                print i,'\t',heads[i]

        files.close()

        return 0

    def printHist(self,filename):

        files=fits.open(filename)

        heads=files[0].header

        if heads['HISTORY']:
            print filename+'\n'
            print heads['HISTORY']
        else:
            print('\n\t --- WARNING: History is not in header ---')

        return 0

    def cleanHead(self,fileName):

        base = fits.open(fileName)
        heads = base[0].header
        
        datas = base[0].data
        datas = np.squeeze(datas)

        if 'NAXIS3' in heads:
            del heads['NAXIS3']
        if 'CRVAL3' in heads:
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

        heads['NAXIS'] = 2

        fits.writeto(fileName,datas,heads,overwrite=True)

        return heads, datas

    def putHead(self,fileName,key,value):

        files=fits.open(fileName)

        datas=files[0].data
        heads=files[0].header

        heads[key] = value

        fits.writeto(fileName,datas,heads,overwrite=True)

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

        add("-hp", "--headPlay",  action="store_true",
                default=True,
                help="call class command")

        add("-prtHd", "--printHead",  action="store_true",
                help="Print header of fits file")

        add("-prtHis", "--printHist",  action="store_true",
                help="Print history of fits file")

        add("-clHd", "--cleanHead",  action="store_true",
                help="Clean header for continuum images and wcs coordinates")

        add("-putHd", "--putHead",  action="store_true",
                help="Put new keyword in header of fits file")

        add('-i', '--input',
            type=str,
            default=False,
            help='''input .fits file''')

        add('-k', '--key',
            type=str,
            default=False,
            help='name of keyword to insert')

        add('-val', '--value',
            type=str,
            default=False,
            help='value of new keyword')

        args = parser.parse_args(argv)

        if args.help:
            print('\n\t************* --- radiobs : headPlay : Help --- **************\n')

            parser.print_help()

            print ("""\nRun the following command:\n
        radiobs -hp -prtHd  -i inputFile.fits 
        radiobs -hp -prtHis -i inputFile.fits 
        radiobs -hp -clHd   -i inputFile.fits 
        radiobs -hp -putHd  -i inputFile.fits  -k <key> -val <value> 

                """)
            print('\n\t************* --- radiobs : headPlay : DONE --- **************\n')

            sys.exit(0)

        elif args.printHead:
            print('\n\t************* ---     radiobs : printHead    --- **************')

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
            self.cleanHead(filename)
            print('\n\t************* ---       Header is clean      --- **************')
            print('\n\t************* --- radiobs : cleanHead : DONE --- **************\n')

        elif args.putHead:
            print('\n\t************* --- radiobs : putHead  --- **************')
            filename = args.input
            self.putHead(filename)
            print('\n\t************* ---   Keword set in header   --- **************')
            print('\n\t************* --- radiobs : putHead : DONE --- **************\n')

        else:
            print ("""\n\t --- ERROR: Please specify input fits file and command ---""")
            print('\n\t************* --- radiobs : headPlay : DONE --- **************\n')
            sys.exit(0)

        print('\n\t************* --- radiobs : headPlay : DONE --- **************\n')



