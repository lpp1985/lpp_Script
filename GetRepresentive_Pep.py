#!/usr/bin/env python
#coding:utf-8
from lpp import *
END = open( sys.argv[2],'w'  )
has = Ddict()


for t,s in fasta_check( open(sys.argv[1],'rU')  ):
	gene = t.split()[0].rsplit("_",1)[0]
	has[gene][len(s)]= t+s
for key in has:
	for data in sorted( has[key] )[::-1]:
		END.write( has[key][ data ]  )
		break
