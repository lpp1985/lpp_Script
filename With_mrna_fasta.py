#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/7/13
from lpp import *
# usage python2.7 RAW_FILE   TARGET_FILE
from lpp import *

from optparse import OptionParser
usage = '''usage: python %prog [options] 

It can automaticly extract fasta with mrna structure!!'''
parser = OptionParser(usage =usage ) 
parser.add_option("-i", "--INPUT", action="store", 
                  dest="input",
                  type='string',
                  help="the input fasta file")

parser.add_option("-m", "--MRNA", action="store", 
                  dest="mrna",
                  type='string',  
                  help="mrna file")

parser.add_option("-o", "--OUTPUT", action="store", 
                  dest="output",
                  type='string',  
                  help="the output file")


(options, args) = parser.parse_args() 


RAW = fasta_check(  open(options.mrna,'rU')  )
FASTA = fasta_check(  open(options.input,'rU')  )
END = open(options.output,'w')
has_hash = {}
for t,s in RAW:
	title = re.search( '>(\d+)',t ).group(1)
	has_hash[  title ] = ''
print( 'The fasta with mrna is %s!!!!'%(  len(  has_hash ) )  )
writ = {}
for t,s in FASTA:
	title = re.search( '>(\d+)',t ).group(1)
	if title in has_hash:
		writ[ title  ] = ''
		END.write( t+s )
print( 'The fasta with output is %s!!!!'%(  len(  writ ) )  )