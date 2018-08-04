#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/7/24
"""

from lpp import * 



if __name__ == '__main__':
	RAW = open(sys.argv[1], 'rU')
	RAW.next()
	END = open(sys.argv[2], 'w')
	for line in RAW:
		line_l = line.split()
		END.write("\t".join( [line_l[3], line_l[2], "similarity",  line_l[4], line_l[5], line_l[1], line_l[6], ".", line_l[8]] )+ '\n')
		
		