#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/7/12
from lpp import *
from optparse import OptionParser

usage = '''usage: python %prog [options] 

It can automaticly cut reads data!!'''
parser = OptionParser(usage =usage ) 
parser.add_option("-i", "--input", action="store", 
                  dest="input",
                  type='string',
                  help="the RAW FASTQ")

parser.add_option("-o", "--output", action="store", 
                  dest="output",
                  type='string',
                  help="the output_FASTQ")

parser.add_option("-s", "--START", action="store", 
                  dest="start",
                  type='int',
                  help= "Start location of reads")
parser.add_option("-e", "--END", action="store", 
                  dest="end",
                  type='int',
                  help="END location of reads ")


(options, args) = parser.parse_args() 




input_fastq = options.input
output_fastq = options.output
start = options.start
end  = options.end

RAW = fastq_check( open( input_fastq , 'rU'  )   )
OUTPUT = open(  output_fastq  ,'w' )
for a,b,c,d in RAW:
	b = b[:-1]
	b_new = b[start:end]
	d = d[:-1]
	d_new = d[ start :  end  ]
	OUTPUT.write( a+b_new+'\n'+c+d_new+'\n'  )
	
	
	