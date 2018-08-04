#!/usr/bin/env python
#coding:utf-8
# Author:  LPP
# Purpose: 
# Created: 2011/11/7
from lpp import *

if __name__=='__main__':
	from optparse import OptionParser 
	usage = '''usage: python2.7 %prog -s sequence_result( 454,or 3730  ) -l gap_location -b blast_result'''
	parser = OptionParser(usage =usage )
	ok_closure = {}
	
	parser.add_option("-q", "--FASTQ", action="store", 
		          dest="fastq", 
		          help="Fasq File")
	
	parser.add_option("-b", "--block_number", action="store", 
		          dest="block", 
	                  type= "int",
		          help="block_number")

	parser.add_option("-o", "--OUTPUT", action="store", 
		          dest="output", 
		          help="prefix of output!!")
	(options, args) = parser.parse_args() 
	seq_file   = options.fastq
	
	block_number  = options.block
	
	output_prefix  = options.output
	
	output_store = {}
	for each_tag in xrange(0, block_number ):
		output_store[ each_tag  ] = open(  '%s.%s.fastq'%(   output_prefix , each_tag  ),'w'  )
	i=0
	for a,b,c,d in fastq_check(  open(   seq_file , 'rU' )  ):
		LOCA = output_store[   i% block_number   ]
		LOCA.write(  a+b+c +d )
		i+=1