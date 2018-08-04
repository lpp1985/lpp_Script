#!/usr/bin/python
from lpp import *
RAW = open(sys.argv[1],'rU')
END = open(sys.argv[2],'w')
for line in RAW:
	line_l = line.split("\t")
	out_l = line_l[:6]
	out_l.append(line_l[1])
	out_l.append(line_l[2])
	out_l.append("0")
	out_l.append("1")
	out_l.append(  str(  int(line_l[2]) -int(line_l[1]) )   )
	out_l.append(  "0"  )
	END.write('\t'.join(   out_l )+'\n'  )
	