#!/usr/bin/env python

import sys, os
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

        print '\n\t************* --- Header is Over --- **************\n'

        print '\n\t************* --- radiobs : printHead : DONE --- **************\n'

    def printHist(self,filename):

        files=fits.open(filename)

        heads=files[0].header

        if heads['HISTORY']:
            print filename+'\n'
            print heads['HISTORY']
        else:
            print '\n\t --- WARNING: History is not in header ---\n'

        print '\n\t************* --- History is Over --- **************\n'
        print '\n\t************* --- radiobs : printHist : DONE --- **************\n'

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

        add('-i', '--input',
            type=str,
            default=False,
            help='''input .fits file''')

        args = parser.parse_args(argv)

        if args.help:
            print '\n\t************* --- radiobs : headPlay : Help --- **************\n'

            parser.print_help()

            print ("""\nRun the following command:\n
        radiobs -hp -prtHd -i inputFile.fits 
        radiobs -hp -prtHis -i inputFile.fits 
                """)
            print '\n\t************* --- radiobs : headPlay : DONE --- **************\n'

            sys.exit(0)

        elif args.prtHd:
            print '\n\t************* --- radiobs : printHead  --- **************\n'

            filename = args.input
            self.printHead(filename)

        elif args.prtHis:
            filename = args.input
            self.printHist(filename)

        else:
            print ("""\n\t --- ERROR: Please specify input fits file ---           """)
            print '\n\t************* --- radiobs : headPlay : DONE --- **************\n'
            sys.exit(0)




