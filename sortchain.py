#!/usr/bin/env python
#coding:utf-8
from lpp import *
RAW =open(sys.argv[1],'rU')
data = Ddict()
for line in RAW:
	if line.startswith("Contig"):
		continue
	line_l = line.split()
	data[line_l[0]][int(line_l[1])]=line

END = open(sys.argv[2],'w')
for ref in data:
	for coord in sorted(data[ref]  ):
		END.write( data[ref][coord]   )
