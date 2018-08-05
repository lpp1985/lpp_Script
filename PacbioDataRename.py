#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/11/4
"""
from lpp import *
all_data  = glob.glob("*.xml")
for f in all_data:
	name = f.split(".")[0]
	result_path = name+'/Analysis_Results/'
	if not os.path.exists(result_path):
		os.makedirs(result_path)
	#print(  "mv %s*.h5 %s "%(name,result_path)    )
	#print(  "mv %s %s "%(f,name)   )
	os.system(  "mv %s*.h5 %s "%(name,result_path)  )
	os.system(  "mv %s %s "%(f,name)  )
