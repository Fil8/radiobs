# Import modules
import os
import sys
import string
import numpy as np


import argparse
from  argparse import ArgumentParser
import textwrap as _textwrap


# get radiobs install directory
RADIOBS_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RADIOBS_DIR = RADIOBS_PATH+'/ants/'
sys.path.append(os.path.join(RADIOBS_PATH, 'ants'))

import pbCorr
import headPlay

import pkg_resources
try:
    __version__ = pkg_resources.require("radiobs")[0].version
except pkg_resources.DistributionNotFound:
    __version__ = "dev"

####################################################################################################

class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width, initial_indent=indent, subsequent_indent=indent) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text

def main (argv):

    for i, arg in enumerate(argv):
        if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

    parser = ArgumentParser(description='radiobs: tools to analyse radio astronomical observations'
                            '|n version {:s} |n install path {:s} |n '
                            'Filippo Maccagni <filippo.maccagni@gmial.com>'.format(__version__,
                                                                               os.path.dirname(__file__)),
                            formatter_class=MultilineFormatter,
                            add_help=False)

    add = parser.add_argument

    add("-h", "--help",  action="store_true",
            help="Print help message and exit")

    add("-v","--version", action='version',
            version='{:s} version {:s}'.format(parser.prog, __version__))

    add('-pb', '--pbCorr',
        action='store_true',
        help= 'primary beam correction')

    add('-hp', '--headPlay',
        action='store_true',
        help= 'tool to play with header of fits file')

    args = parser.parse_args(argv)

    if args.help and len(argv) ==1 :
        print '\n\t************* --- radiobs : Help --- **************\n'

        print ('\t\t  ... called for help ...\n')
        parser.print_help()

        print ("""\nRun a command. This can be:\n
radiobs\t\t(all tools)
radiobs -pb\t correction for the primary beam
            """)
        print '\n\t************* --- radiobs : DONE --- **************\n'

        sys.exit(0)

    elif args.pbCorr:
        
        print ('\n\t************* --- radiobs : pbCorr --- **************\n')
        pb = pbCorr.pbcorr()
        pb.main(argv)

    elif args.headPlay:
        
        print ('\n\t************* --- radiobs : headPlay --- **************\n')
        hp = headPlay.headplay()
        hp.main(argv)

    else:
        print ('\n\t ... you have not entered an available class function ... \n')
        print ('\t************* --- radiobs : ERROR --- **************\n')
        sys.exit(0)