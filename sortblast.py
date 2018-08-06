#!/usr/bin/env python
#coding:utf-8
from lpp import *
RAW = open(sys.argv[1],'rU')
END = open( sys.argv[2],'w' )
has = {}
already= []
for line in RAW:
	line_l = line.split("\t")
	if len(has)==0:
		has[line_l[0]]=""
	if line_l[0] not in has:
		has[line_l[0]]=""
		already = sorted(already,key=lambda x: int(x.split("\t")[6]) )
		END.write(''.join(already))
		already = []
	already.append(line)
already = sorted(already,key=lambda x: int(x.split("\t")[6]) )
END.write(''.join(already))
