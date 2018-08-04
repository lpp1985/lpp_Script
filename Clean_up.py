#!/usr/bin/env python
#coding:gb2312
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/4/6
from lpp import *
import sys, os, zipfile,shutil,glob,os
from optparse import OptionParser 

usage = '''usage: python2.7 %prog -i input_path -o output_path'''
parser = OptionParser(usage =usage ) 
parser.add_option("-i", "--INPUT", action="store", 
                  dest="input",
                  default = './', 
                  help="Input Path")
parser.add_option("-a", "--append", action="store", 
                  dest="append", 
                  help="appendix")
(options, args) = parser.parse_args() 
append = options.append
intputpath = options.input
def check_path( path ):
	if not os.path.exists(path):
		os.mkdir( path )

for each_f in glob.glob( '*.%s'%( append  ) ):
	title = each_f.split( '.' )[0]
	check_path( title )
	for e_f in glob.glob(  '%s.*'%( title ) ):
		shutil.move( e_f,title +'/' )