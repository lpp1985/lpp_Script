#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/1/11
"""
import sys,copy
from optparse import OptionParser


if __name__ == '__main__':
	usage='''usage: python %prog [options]

        Transform RepeatMasker gff result'''
	parser = OptionParser(usage =usage )
	parser.add_option(
		"-i", "--Input", action="store",
		dest="out",
		type='string',
		help="RepeatMasker out file")


	parser.add_option(
		"-o", "--out", action="store",
		dest="output",
		type='string',
		help="Output Gff File")     
	(options, args) = parser.parse_args()
	data_hash = {}

	REPOUT = open( options.out,'rU'  )
	RESULT = open( options.output,'w')

	for line in REPOUT:
		if line.startswith("#"):
			RESULT.write(line)
			continue
		line_l = line.split("\t")
		new_l = copy.deepcopy(line_l)
		new_l[1] = line_l[3]
		new_l[2] = line_l[4]
		new_l[4] = line_l[2]
		new_l[3] = line_l[1]
		RESULT.write("\t".join(new_l) )
