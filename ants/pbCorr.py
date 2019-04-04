#!/usr/bin/env python

import os, sys, string
import numpy as np
from astropy.io import fits

import argparse
from  argparse import ArgumentParser
import textwrap as _textwrap

import pbCorr

class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width, initial_indent=indent, subsequent_indent=indent) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text

class pbcorr:

    def __init__(self):

        self.rootdir = os.getcwd()+'/'


    def BeamCorrect(self,fileName,datas,hdr,telescope):

        try :
            hdr['CRVAL3']
        except KeyError:
            obs_freq = float(hdr['FREQ'])   
        else:
            obs_freq = float(hdr['CRVAL3'])

        if 'RESTFREQ' in hdr:
            hdr['FREQ'] = hdr['RESTFREQ']

        if telescope == 'MeerKAT':
            ant = 13.5
        elif telescope == 'VLA':
            ant = 25.
        elif telescope == 'ACA':
            ant = 7.

        pb_fwhm = 1.02*(2.99792458E8)/obs_freq/ant/np.pi*180.
        pb_fwhm_pix = pb_fwhm/hdr['CDELT2']
        x, y = np.meshgrid(np.linspace(-hdr['NAXIS1']/2.,hdr['NAXIS1']/2.,hdr['NAXIS1']), 
                           np.linspace(-hdr['NAXIS2']/2.,hdr['NAXIS2']/2.,hdr['NAXIS2']))
        d = np.sqrt(x*x+y*y)
        sigma, mu = pb_fwhm_pix/2.35482, 0.0
        gauss = np.exp(-( (d-mu)**2 / ( 2.0 * sigma**2 ) ) )
        
        pbcorr = np.divide(datas,gauss)
        outfile=string.split(fileName,'.fits')[0]
        #if len(string.split(outfile,'/')) >1:
    #       outfile=string.split(outfile,'/')[1]
        
        out_pbcorr = outfile+'_pbcorr.fits'
        out_beam = outfile+'_gauss_beam.fits'

        fits.writeto(out_beam,gauss,hdr,overwrite=True) 
        fits.writeto(out_pbcorr,pbcorr,hdr,overwrite=True)


        return out_pbcorr


    def main (self,argv):

        for i, arg in enumerate(argv):
            if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

        parser = ArgumentParser(description='radiobs: pbCorr : primary beam correction',
                                formatter_class=MultilineFormatter,
                                add_help=False)

        add = parser.add_argument

        add("-h", "--help",  action="store_true",
                help="Print help message and exit")

        add('-i', '--input',
            type=str,
            default=False,
            help='''input .fits file''')

        add('-tel', '--telescope',
            type=str,
            default=False,
            help='select telescope: MeerKAT, VLA, ACA')

        args, unknown = parser.parse_known_args()

        if args.help:
            print '\n\t************* --- radiobs : pbCorr : Help --- **************\n'

            parser.print_help()

            print ("""\nRun the following command:\n
        radiobs -pb -i inputFile.fits -tel <MeerKAT,VLA,ACA>
                """)
            print '\n\t************* --- radiobs : pbCorr : DONE --- **************\n'

            sys.exit(0)

        elif args.input:
            
            if args.telescope == False:

                print ("""\n\t --- ERROR: Please specify telescope : MeerKAT,VLA,ACA> ---
                    """)
                print '\n\t************* --- radiobs : pbCorr : DONE --- **************\n'
                sys.exit(0)

            pb = pbCorr.pbcorr()
            fileName=args.input
            telescope=args.telescope
            ff=fits.open(fileName)
            dats=ff[0].data
            heads=ff[0].header
            dats=np.squeeze(dats)
            pb.BeamCorrect(fileName,dats,heads,telescope)  

            print '\n\t************* --- Primary Beam Correction --- **************\n'

            print '\n\t************* --- radiobs : pbCorr : DONE --- **************\n'

        else:

            print ("""\n\t --- ERROR: Please specify input fits file ---
                    """)
            print '\n\t************* --- radiobs : pbCorr : DONE --- **************\n'
            sys.exit(0)
