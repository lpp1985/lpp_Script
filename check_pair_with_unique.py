#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/4
from lpp import *
from optparse import OptionParser 
usage = '''usage: python2.7 %prog -1 read1 -2 read2 -o output'''
parser = OptionParser(usage =usage ) 
parser.add_option("-1", "--Reads1", action="store", 
                  dest="read1", 
                  help="Read1 location") 
parser.add_option("-2", "--Reads2", action="store", 
                  dest="read2", 
                  help="Read2 location")

parser.add_option("-o", "--Output", action="store", 
                  dest="output", 
                  help="output")

parser.add_option("-u", "--Unpaired", action="store", 
                  dest="unpair", 
                  help="unpair")
(options, args) = parser.parse_args() 
read1 = options.read1
read2 = options.read2
output = options.output
R1 = fastq_check( open( read1,'rU'  )   )

UNPAIR = open( options.unpair,'w'  )

r1_have = {}
for t,s,w,q in R1:
	title = t[:-2]
	r1_have[ title  ] = t+s+w+q
R2 = fastq_check( open( read2,'rU'  )   )
R1_end = open( output+'.read1.fastq','w'  )
R2_end = open( output+'.read2.fastq','w'  )
all_have = {}
r2_have = {}
for t,s, w,q in R2:
	title = t[:-2]
	r2_have[ title ] = ''
	if title in r1_have:
		R2_end.write( t+s+w+q )
		R1_end.write( r1_have[title]  )
		all_have[ title ] = ''

	else:
		UNPAIR.write( title+'\tR2_have\n' )
		
for each_key in r1_have:
	if each_key not in r2_have:
		UNPAIR.write( title+'\tR1_have\n' )