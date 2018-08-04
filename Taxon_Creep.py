#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/6/9
"""
from lpp import *
RAW = open(sys.argv[1],'rU')
all_spieces = Ddict()
for line in RAW:
    line_l = line.split("\t")
    if line_l[2]!=line_l[0]:
        all_spieces[line_l[2]][line_l[0]]=""
        
def get_taxon_seed2( taxon_number ):
    all_end={}
    def creeper(taxon_number):
	all_end[ taxon_number ]=''
        if taxon_number in all_spieces:
            for key1 in all_spieces[ taxon_number ]:
                creeper(key1)
    creeper(taxon_number)
    return all_end

END = open(sys.argv[3],'a')
data = sys.argv[2]
all_need = get_taxon_seed2(data)

for key in all_need:
        END.write(key+"\t"+data+'\n')
