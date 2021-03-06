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

argv = [a for a in sys.argv[1:]]

for i, arg in enumerate(argv):
    if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

parser = ArgumentParser(description='radiobs: PbCorr primary beam correction',
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

args = parser.parse_args(argv)

if args.help:
    print '\n\t************* --- radiobs : PbCorr : Help --- **************\n'

    parser.print_help()

    print ("""\nRun the following command:\n
pbCorr -i inputFile.fits -tel <MeerKAT,VLA,ACA>
        """)
    print '\n\t************* --- radiobs : PbCorr : DONE --- **************\n'

    sys.exit(0)

elif args.input:
    
    if args.telescope == False:

        print ("""\n\t --- ERROR: Please specify telescope : MeerKAT,VLA,ACA> ---
            """)
        print '\n\t************* --- radiobs : PbCorr : DONE --- **************\n'
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

    print '\n\t************* --- radiobs : PbCorr : DONE --- **************\n'

