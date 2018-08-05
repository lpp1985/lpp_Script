#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/9/26
"""

from lpp import *
RAW = open(sys.argv[1],'rU')
END = open(sys.argv[1]+'1','w')
for line in RAW:
	line_l = line.split("\t")
	line_l[-1]=line_l[-1].replace(".",'_')
	END.write("\t".join(line_l))
