#!/usr/bin/python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/6/3
from lpp import *
all_data = Ddict()
all_file = sys.argv[1:-1]
RESULT = open(sys.argv[-1],'w')
all_status = {}
for e_f in all_file:
	RAW = open(e_f,'rU')
	title = RAW.next()
	for line in RAW:
		line_l = line.strip().rsplit("\t",2)
		all_data[ line_l[0] ][ line_l[-1]   ]=line_l[-2]
		all_status[ line_l[-1] ] = ""
RESULT.write(title)
for key in all_data:
	for status in all_status:
		if status not in all_data[key]:
			RESULT.write(key+'\t0\t'+status+'\n')
		else:
			RESULT.write( key+'\t'+all_data[key][status]+'\t'+status+'\n'   )
