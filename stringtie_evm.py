#!/usr/bin/env python
#coding:utf-8
"""
  Author:   -->
  Purpose: 
  Created: 2015/10/24
"""
from lpp import *
RAW = open(sys.argv[1],'rU')
END = open(sys.argv[2],'w')
data_hash = {}
data_all = RAW.read()
data_block = re.split(".+\ttranscript\t.+\n",data_all)
i=0
for each_block  in data_block[1:]:
    i+=1
    total_length = 0
    each_block = each_block.strip()
    if not each_block:
        continue
    for line in each_block.split("\n"):
        line_l = line.split("\t")
        line_l[0] = line_l[0].split("|")[0]
        line_l[2] = "EST_match"
        new_data = line_l[:2]    
        length = int(line_l[4]) - int(line_l[3])+1
        total_length+=length
        gene_id = re.search("gene_id \"(\S+)\"\;",line).group(1)
        transcript_id=re.search("transcript_id \"(\S+)\"\;",line).group(1)
    start =1    
    
    for line in each_block.split("\n"):
        
        line_l = line.split("\t")
        line_l[0] = line_l[0].split("|")[0]
        line_l[2] = "EST_match"
        line_l[4] = str(int(line_l[4])+1)
        align_length = int(line_l[4]) -int(line_l[3])
        line_l[5] = "%.2f"%(100.0*align_length/total_length)
        line_l = line_l[:-1]

        line_l.append ("ID=Align_%s;Target=EST_%s %s %s"%(i,i,start,start+align_length) )
        start = start+align_length+1
        END.write('\t'.join(line_l)+'\n')