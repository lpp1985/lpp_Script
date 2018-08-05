#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2016/5/17
"""

from lpp import *
RAW = pd.read_table(sys.argv[1])
all_counts_name  = RAW.columns[6:]
for key in all_counts_name:
    new_data = pd.DataFrame( RAW,columns=[ "Geneid",key ]  )
    new_data.to_csv(key+'.count',header=False,index=False,sep="\t")
    
    