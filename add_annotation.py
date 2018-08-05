#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/6/19
"""
from lpp import *
RAW = fasta_check(open(sys.argv[1],'rU') )
anno_hash = {}
for t,s in RAW:
	t = t[1:-1]
	i_d,name= t.split(' ',1)
	anno_hash[i_d] = name
#print(anno_hash)
	
END = open("SNP_anno.list",'w')
for line in open(sys.argv[2],'rU'):
	line_l = line.split("\t")
	if line_l[0] in anno_hash:
		END.write(line[:-1]+'\t'+anno_hash[line_l[0]]+'\n')
	else:
		END.write(line)
