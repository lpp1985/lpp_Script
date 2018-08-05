#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/11/3
"""

from lpp import *
all_name = File_dict(open(sys.argv[1], 'rU')).read(1, 1)
all_data = set(all_name)

if __name__ == '__main__':
	END = open("SingleOrtholog.list", 'w')
	for line in open(sys.argv[2], 'rU'):
		line_l = line.split(": ")
		group_name = line_l[0]
		all_gene = line_l[1].strip().split()
		if len(all_gene) == len(all_data):
			all_org = set(re.findall("(\S+)\|", line_l[1]))
			if all_org == all_data:
				END.write(line_l[0] + '\t' + line_l[1])
				
			
		