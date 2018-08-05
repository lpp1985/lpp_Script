#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2018/1/17
"""

from lpp import * 



if __name__ == '__main__':
	SEQ = fasta_check( open(sys.argv[1]))

	for t, s in SEQ:
		all_seq = re.sub("\s+", "", s)
	DATA = open(sys.argv[2])
	END = open(sys.argv[3], 'w')
	num = 0
	for line in DATA:
		if ".." in line:
			data = line.split()[2]
			start , end = data.split("..")
			start = int(start)
			end = int(end)
			num += 1
			END.write('>Alien%s\n' % (num))
			END.write(all_seq[start:end] + '\n')
		