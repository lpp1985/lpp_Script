#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/11/24
"""

from lpp import * 
from optparse import OptionParser
from termcolor import colored
import os


if __name__ == '__main__':
	
	usage = '''usage: python %prog [options] 
	
	It can automaticly run BWA!!'''
	parser = OptionParser(usage =usage ) 
	
	parser.add_option("-i", "--Input", action="store", 
		              dest= "Input",
		              type='string',
		              help="the input seq")
	
	parser.add_option("-o", "--Output", action="store", 
		              dest="output",
		              type='string',
		              help="output")
	
	
	(options, args) = parser.parse_args()
	
	Input = options.Input
	output = options.output
	all_result = os.popen( "tmhmm %s" % (Input))
	data = all_result.read()
	RESULT = open( output, 'w')
	data_l = re.split( "\n(?=\#[^\n]+Length\:)", data)
	RESULT.write("Name\tTMHMM_Domain\tTMHMMStart\tTMHMMStop\n")
	for block in data_l:
		number = re.search("Number of predicted TMHs\:\s+(\d+)", block).group(1)
		if number != 0:
			line_all = block.split("\n")
			for line in line_all:
				if line and  not line.startswith ("#"):
					line_l = line.split("\t")
					line_l.pop(1)
					RESULT.write("\t".join(line_l) + '\n')
	os.system("rm TMHMM_*/  -rf ")
