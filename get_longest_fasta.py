#!/usr/bin/env python
#coding:utf-8
"""
  Author:   -->
  Purpose: 
  Created: 2015/12/7
"""
from lpp import *
from optparse import OptionParser

if __name__=="__main__":
	usage = '''usage: python2.7 %prog'''
	parser = OptionParser(usage =usage ) 
	parser.add_option("-i", "--INPUT", action="store", 
	                  dest="input", 
	                  help="input file")

	parser.add_option("-o", "--end", action="store", 
	                  dest="output", 
	                  help="output fasta file")

	parser.add_option("-l", "--length", action="store", 
	                  dest="length",
	                  type="int",
	                  help="evalue cutoff")
	(options, args) = parser.parse_args()
	RAW = fasta_check( open(options.input,'rU')  )
	END = open(options.output, mode='w')
	length_hash = Ddict()
	for t,s in RAW:
		length = len(re.sub("\s+","",s))
		length_hash[length][t+s]=""
	all_length = 0
	for length in sorted(length_hash)[::-1]:
		for key in length_hash[length]:
			all_length+=length
			
			if all_length<options.length:
				END.write(key)
			else:
				sys.exit()
	
