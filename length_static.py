#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/11/4
"""

from lpp import *
RAW = fasta_check(open(sys.argv[1],'rU'))
for t,s in RAW:
	s = re.sub("\s+",'',s)
	print(len(s))
