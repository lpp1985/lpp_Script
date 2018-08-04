#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: Trans 454 to fastq
# Created: 2011/8/8
from lpp import *
from optparse import OptionParser
usage = '''usage: python %prog [options] 

It can automaticly transfer 454 to fastq!!'''
parser = OptionParser(usage =usage ) 
parser.add_option("-s", "--FASTA", action="store", 
                  dest="fasta",
                  type='string',
                  help="the fasta file name")

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


parser.add_option("-t", "--Other", action="store", 
                  dest="left",
                  type='string',  
                  help="THE output of left  fastq file")

(options, args) = parser.parse_args() 

FASTA = fasta_check(  open( options.fasta,'rU'  )   )

LEFT = open( options.left,'w'   )
all_name ={}
RAW =  open( options.list,'rU'   )
for line in RAW:
    line_l = line[:-1].split('\t')
    all_name[  line_l[options.number-1].split()[0]] = ''
END = open( options.output ,'w' )
for t,s in FASTA:
	title = t[1:-1].split()[0]
	if title in all_name:
		END.write(  t+s  )
	else:
		LEFT.write( t+s  )
		
		
		




