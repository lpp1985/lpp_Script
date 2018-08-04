#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/1/11
"""
from lpp import *
END = open(sys.argv[2],'w')
RAW = open(sys.argv[1] ,'rU')
RAW.next()
END.write('Category\tAnnotation\tNumber\n')
cog_mapping = File_Ddict( open(sys.argv[1],'rU') ).read(7,1)
cog_detail = File_dict( RAW  ).read(7,8)
for key ,val in cog_detail.items():
	#print(cog_mapping[key])
	END.write(key+'\t'+val +'\t%s\n'%(   len(cog_mapping[key])))
