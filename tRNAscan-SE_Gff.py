#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/5/21
"""

from lpp import *
RAW = open(sys.argv[1])
RAW.next()
RAW.next()
RAW.next()
END = open(sys.argv[2],"w")
gene = 1
for line in RAW:
	line_l = line.split()
	contig = line_l[0]
	begin = line_l[2]
	end = line_l[3]
	product = line_l[4]
	if int(begin)> int(end):
		frame = '-'
		begin,end = end,begin
	else:
		frame="+"
		
	END.write("%s\ttRNAscan-SE\tgene\t%s\t%s\t.\t%s\t.\tID=tRNA%s\n"%(contig,begin,end,frame,gene))
	END.write("%s\ttRNAscan-SE\ttRNA\t%s\t%s\t.\t%s\t.\tID=tRNA%s;Parent=tRNA%s;product=tRNA-%s\n"%(contig,begin,end,frame,gene,gene,product))
	END.write("%s\ttRNAscan-SE\texon\t%s\t%s\t.\t%s\t.\tID=tRNA%s;Parent=tRNA%s;product=tRNA-%s\n"%(contig,begin,end,frame,gene,gene,product))
	gene+=1