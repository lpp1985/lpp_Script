#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/11/13
"""

from lpp import *
all_member = File_Ddict(open(sys.argv[1],'rU')).read(2,1)
ALL_SEQ = fasta_check(open(sys.argv[2],'rU'))

if __name__ == '__main__':
	BAC_SEQ = open("All_Bacteria.fasta",'w')
	for t,s in ALL_SEQ:
		if t[1:-1] in all_member:
			BAC_SEQ.write(t+s)