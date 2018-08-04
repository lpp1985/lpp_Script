#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/4/24
"""

from lpp import *
RAW = pd.read_table(sys.argv[1])
EXCEL= pd.ExcelWriter(sys.argv[2], engine='xlsxwriter')
RAW.to_excel( EXCEL,"Sheet1",index=False    )