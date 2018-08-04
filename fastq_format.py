#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/8/21
"""

from lpp import *
RAW = fastq_check(open(sys.argv[1],'rU'))
END = open(sys.argv[2],'w')
for a,b,c,d in RAW:
    END.write(a+b+'+'+a[1:]+d)
