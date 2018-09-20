#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2018/9/20
"""

from lpp import *
import sys



if __name__ == '__main__':
    all_data = Ddict()
    for line in sys.stdin:
        line_l = line.strip().split("\t")
	if len(line_l)<5:
		continue
	if line_l[4] not in ["-","+"]:
		continue
        all_data[line_l[0]][int(line_l[-1])]=line_l
    for seq in all_data:
        need = sorted(all_data[seq].keys())[-1]
        all_data[seq][need].append( "%.2f"%(100.0*( int(all_data[seq][need][3])-  int(all_data[seq][need][2])  ) / int(all_data[seq][need][1]))  )
	print("\t".join(all_data[seq][need]))

