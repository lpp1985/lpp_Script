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


parser.add_option("-i", "--INPUT", action="store", 
                  dest="input_data",
                  type='string',  
                  help="the  file name to input")


parser.add_option("-o", "--OUTPUT", action="store", 
                  dest="output",
                  type='string',  
                  help="THE output fastq file")



(options, args) = parser.parse_args() 

RAW = fastq_check(   open(  options.input_data,'rU'  )  )
END = open( options.output,'w' )
for a,b,c,d in RAW:
	b = complement(  b[:-1] )+'\n'
	d = d[:-1][::-1]+'\n'
	END.write( a+b+c+d  )
		
		
		




