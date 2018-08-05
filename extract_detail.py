#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: Trans 454 to fastq
# Created: 2011/8/8
from lpp import *
from optparse import OptionParser
usage = '''usage: python %prog [options] 

It can automaticly get information you want!!'''
parser = OptionParser(usage =usage ) 
parser.add_option("-s", "--Input", action="store", 
                  dest="input_data",
                  type='string',
                  help="the input  file name")

parser.add_option("-l", "--LIST", action="store", 
                  dest="list",
                  type='string',  
                  help="the list  file name to fetch")

parser.add_option("-n", "--number", action="store", 
                  dest="number",
                  type='int',  
                  help="the number of column to extract [ from 1 to start   ]")

parser.add_option("-o", "--OUTPUT", action="store", 
                  dest="output",
                  type='string',  
                  help="THE output  file")
parser.add_option("-a", "--ANOTHER", action="store", 
                  dest="another",
                  type='int',  
                  help="The column number in the data file(from 1 start)")

(options, args) = parser.parse_args() 

RAW = open( options.input_data,'rU'  )   
another = options.another

all_name = File_Ddict(  open( options.list,'rU'   )   ).read( options.number , options.number   )


END = open( options.output ,'w' )
for line in RAW:
	title = line.strip().split('\t')[another-1]
	if title in all_name:
		END.write(  line  )

		
		
		




