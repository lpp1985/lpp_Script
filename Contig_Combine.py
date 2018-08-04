#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/8/25
"""

from lpp import * 



if __name__ == '__main__':
	RAW = fasta_check(open(sys.argv[1], 'rU'))
	END = open(sys.argv[2], 'w')
	seq = ""
	for t, s in RAW:
		seq += "N" * 20 + s
	seq = re.sub("\s+", "", seq)
	seq = re.sub("N{20,}", "N"*20, seq)
	seq = re.sub("^N+","",seq)
	seq = re.sub("N+$","",seq)
	seq = re.sub("(\w{60})", "\\1\n", seq)
	END.write(">Scaffold\n%s\n" % (seq))
		
