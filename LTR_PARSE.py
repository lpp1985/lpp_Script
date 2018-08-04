#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/8/18
"""

from lpp import * 
from optparse import OptionParser


if __name__ == '__main__':
	usage = '''usage: python2.7 %prog [options] Kmer
	
	
	
	
		    Kmer is a list of K value you want,e.g  [ 1, 2, 3, 4 ]'''
	parser = OptionParser(usage =usage )



	parser.add_option("-i", "--Input", action="store",
                          dest="inp",

	                      help="input data")
	parser.add_option("-o", "--output", action="store",
                          dest="out",
	                      
	                      help="output")
	(options, args) = parser.parse_args()
	ltr_raw = options.inp
	OUTPUT = open( options.out, 'w')
	RAW = fasta_check( open(ltr_raw) )
	for t, s in RAW:
		for line in s.split("\n"):
			if line.startswith("["):
				line = re.split("\]\s+",line)[1]
				line_l = line.split()
				start, end = line_l[1].split("-")
				OUTPUT.write( "\t".join([ line_l[0], "LTR-FInder", "Repeat", start, end,  ".", "+", ".", "LTR_Repeat" ] )+ '\n' )
				
		
		
