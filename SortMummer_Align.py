#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2018/2/27
"""

import sys


if __name__ == '__main__':
	data = sys.stdin
	all_data = []
	print(data.next()), 
	for line in data:
		all_data.append(line)
	all_data = sorted(all_data, key = lambda x: int(x.split("\t")[4]) )
	all_data = sorted(all_data, key = lambda x: x.split("\t")[1]) 
	for key in all_data:
		
		print(key), 
