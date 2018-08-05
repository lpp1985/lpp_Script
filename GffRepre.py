#!/usr/bin/env python
#coding:utf-8
from lpp import *
END = open( sys.argv[2],'w'  )
has = {}
for t,s in fasta_check( open(sys.argv[1],'rU')  ):
	gene = t.split()[1]
	if gene not in has:
		has[gene]=''
		END.write(t+s)
