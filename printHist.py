#!/usr/bin/env python

from astropy.io import fits
import sys

import argparse
from  argparse import ArgumentParser
import textwrap as _textwrap

argv = [a for a in sys.argv[1:]]

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

parser = ArgumentParser(description='radiobs: prtHist print history of fits file',
                        formatter_class=MultilineFormatter,
                        add_help=False)

add = parser.add_argument

add("-h", "--help",  action="store_true",
        help="Print help message and exit")

add('-i', '--input',
    type=str,
    default=False,
    help='''input .fits file''')


args = parser.parse_args(argv)

if args.help:
    print '\n\t************* --- radiobs : prtHist : Help --- **************\n'

    parser.print_help()

    print ("""\nRun the following command:\n
prtHist -i inputFile.fits 
        """)
    print '\n\t************* --- radiobs : prtHist : DONE --- **************\n'

    sys.exit(0)

elif args.input:

    filename = args.input
    files=fits.open(filename)

    heads=files[0].header

    if heads['HISTORY']:
        print filename+'\n'
        print heads['HISTORY']
    else:
        print '\n\t --- WARNING: History is not in header ---\n'

else:

        print ("""\n\t --- ERROR: Please specify input fits file ---
                    """)
        print '\n\t************* --- radiobs : prtHist : DONE --- **************\n'
        sys.exit(0)

