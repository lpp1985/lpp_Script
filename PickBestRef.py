#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2019/6/17
"""
from lpp import *


if __name__ == '__main__':
    reads_number = Ddict()
    final_reads = {}
    read_category = {}
    all_seq = {}
    END_BEST = open("Best1.fa",'w')
    for t,s in fasta_check( open(sys.argv[1],'rU') ):
        cat = re.search(  "\s+\((\S+)\)",t ).group(1)
        name = t.split()[0][1:]
        read_category[name]= cat
        all_seq[ t[1:].split()[0] ]=s
        
    for line in open(sys.argv[2],'rU'):
        line_l = line.strip().split("\t")
        reads_number[  read_category[ line_l[0]  ]   ][line_l[0]]=int(line_l[1])
        
        
    for key in reads_number:
        for v in sorted( reads_number[key ],key=lambda x: reads_number[key ][x]  )[::-1]:
            END_BEST.write('>'+v+"_%s"%(key)+"\n"+all_seq[v])
            break
