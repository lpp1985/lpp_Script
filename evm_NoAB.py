#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/11/21
"""

from lpp import *
RAW = open(sys.argv[1],'rU')
END = open(sys.argv[1]+'1','w')

DATA = RAW.read().split("\n\n")
title,data = DATA[0].split("\n",1) 

END.write(title+'\n')
DATA[0] = data
for data in DATA:
    if "StringTie" not in data and "GeneWise" not in data:
        continue
    END.write(data+'\n\n')

    
