#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/3/7
"""
from lpp import *
RAW = fasta_check(open(sys.argv[1],'rU'))
END = open(sys.argv[3],'w')
need_data = File_dict(open(sys.argv[2],'rU')).read(1,1)
for t,s in RAW:
	gi_all = re.findall("gi\|(\d+)", t)
	for i in gi_all:
		if i in need_data:
			END.write(t+s)
			break
		
