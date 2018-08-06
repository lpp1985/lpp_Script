#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2018/7/25
"""

import unittest
from lpp import *
from networkx import DiGraph
from Dependency import *
if __name__ == '__main__':
    RAW = open(sys.argv[1],'rU')
    END = open(sys.argv[2],'w')
    has = {}
    seq = ""
    number = 1
    network = Contig_Graph()
    for line in RAW:
        line_l = line.strip().split()
            
        if line_l[0] not in has:
            has[line_l[0]] = ""
	    start = line_l[-3]+line_l[-2]
                
        else:
            end =  line_l[-3]+line_l[-2]
            network.add_bi_edge(start,end)
            start = end 
    for start ,end in network.edges():
        END.write(start+'\t'+end+'\n')
    #END.write('>scaffold%s\n'%(number)+seq+'\n')
            
