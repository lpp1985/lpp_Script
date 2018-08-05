#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/3/30
"""
from lpp import *
RAW= fasta_check(open(sys.argv[1]))
length_data = Ddict()
for t,s in RAW:
	s =re.sub("\s+","",s).upper()
	length = len(s)
	length_data[length][t] = ""
END = open("length_distribution.tsv","w")
data = 0
for length in reversed(sorted(length_data)):
	data+= length*len(length_data[length])
	END.write("%s\t%s\n"%( length,data   )   )