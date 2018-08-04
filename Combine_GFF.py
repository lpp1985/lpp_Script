#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/1/16
"""

from lpp import *



if __name__ == '__main__':
	result = Ddict()
	END = open(sys.argv[-1],'w')
	for e_f in sys.argv[1:-1]:
		RAW=open(e_f,'rU')
		for line in RAW:
			if "Parent=" not in line:
				line_l = line.split("\t")

				ID = re.search("\tID\=([^\;]+)",line).group(1)
				chro = line_l[0]
				coor = int(line_l[3])
				result[chro][coor][ID] = line
			else:
				result[chro][coor][ID] += line
	for key in sorted(result):
		for coor in sorted(result[key]):
			for data in sorted( result[key][coor] ):
				END.write(result[key][coor][data])