#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/3/13
"""
from lpp import *
from optparse import OptionParser
usage = "python2.7 %prog [options]"
parser = OptionParser(usage =usage )
parser.add_option("-i", "--XLS", action="store",
                  dest="xls",

                  help="IS_Annotation")

parser.add_option("-o", "--out", action="store",
                  dest="output",

                  help="output")




if __name__ == '__main__':
	(options, args) = parser.parse_args()
	DATA = open(   options.xls)
	stat =Ddict()
	DATA.next()
	for line in DATA:
		line_l = line.split("\t")
		if line_l[8]  not in stat or  line_l[4] not in stat[line_l[8]]:
			
			stat[ line_l[8] ][  line_l[4]]=[int(  line_l[6] )]
		else:
			stat[ line_l[8] ][  line_l[4]].append(int(  line_l[6] ))
	
	output = options.output
	OUTPUT = open(  output,'w' )
	OUTPUT.write(   "Family\tGroup\tNumber\tAverage Length\n"  )
	for fam in stat:
		for group in stat[fam]:
			number = len( stat[fam][group]  )
			av_length = sum( stat[fam][group]    )/number
			OUTPUT.write(fam+'\t'+group+'\t'+str(   number  )+'\t'+str(av_length)+'\n'     )
	