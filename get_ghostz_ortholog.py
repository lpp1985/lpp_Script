#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/6/29
"""

from lpp import *
RAW = open(sys.argv[1],'rU')
all_ortholog = File_dict(open(sys.argv[2],'rU')).read(1,1)
HAS = open(sys.argv[3],'w')
OTHER = open(sys.argv[4],'w')
for line in RAW:
    name = line.split("\t")[0].split()[0]
    if name in all_ortholog:
        HAS.write(line)
    else:
        OTHER.write(line)