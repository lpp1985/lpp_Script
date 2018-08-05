#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2016/3/16
"""

import pandas as pd
import sys
raw = pd.read_table(sys.argv[1])
new_frame = pd.DataFrame(raw,columns= sys.argv[2:])
new_frame.to_csv("need.tsv",index=False,sep="\t")