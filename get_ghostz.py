#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/6/27
"""
from lpp import *
RAW = File_dict(open(sys.argv[1],'rU')).read(1,1)
END = open(sys.argv[3],'w')
FROM = open(sys.argv[2],'rU')
for line in FROM:
    line_l = line.split()
    if line_l[0] in RAW:
        END.write(line)
