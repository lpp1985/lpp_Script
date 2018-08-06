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

if __name__ == '__main__':
    RAW = open(sys.argv[1],'rU')
    END = open("ScaffoldEdges.tsv",'w')
    has = {}
    seq = ""
    number = 1
    network = DiGraph()
    for line in RAW:
        line_l = line.strip().split()
            
        if line_l[0] not in has:
            has[line_l[0]] = ""
	    start = line_l[7]+line_l[8]
                
        else:
            end =  line_l[7]+line_l[8]
            network.add_edge(start,end)
            start = end 
    for start ,end in network.edges():
        END.write(start+'\t'+end+'\n')
    #END.write('>scaffold%s\n'%(number)+seq+'\n')
            
