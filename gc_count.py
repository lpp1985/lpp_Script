#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/11/20
"""
from lpp import *
def GC_Count(data):
	FASTA = fasta_check(open(data,'rU'))
	total = 0
	gc =0

	for t,s in FASTA:
		total+= len(re.findall("((?:G|C|A|T))", s))
		gc += len(re.findall("((?:G|C))",s))
		
	return float(gc)/total
print(GC_Count(sys.argv[1]))