#!/usr/bin/python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: Trans 454 to fastq
# Created: 2011/8/8
from lpp import *
from optparse import OptionParser
usage = '''usage: python %prog [options] 

It can automaticly transfer 454 to fastq!!'''
parser = OptionParser(usage =usage ) 
parser.add_option("-m", "--m4", action="store", 
                  dest="m4",
                  type='string',
                  help="the m4 file name")

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
                  help="THE output fastq file")

(options, args) = parser.parse_args() 

M4 =  open( options.m4,'rU'  ) 


all_name = File_Ddict(  open( options.list,'rU'   )   ).read( options.number , options.number   )
all_have = {}
END = open( options.output ,'w' )
for line in M4:
	line_l = line.split("\t")
	title = line_l[0].rstrip().split()[0]
	if title in all_name:
		all_have [ title  ] = ''
		END.write(  line )
for key in all_name :
	if key not in all_have:
		print( key )





