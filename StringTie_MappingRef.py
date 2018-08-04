#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2016/5/16
"""

from lpp import *
import re
RAW = open(sys.argv[1],'rU')
mapping = {}
MAPPING_ALL= open( sys.argv[3],'w')
total_mapping = {}
for line in RAW:
    line_l = line.split("\t")
    if "\ttranscript\t" in line:
        all_ids = re.findall("\S+ \"(\S+)\"", line_l[-1])
        
        gene_name = re.search(  "gene_name \"(\S+)\"",line )
        gene_id = re.search(  "gene_id \"(\S+)\"",line ).group(1)
        raw_id = re.search(  "gene_id \"(\S+)\"",line ).group(1)
            
        ref_gene_id = re.search(  "ref_gene_id \"(\S+)\"",line )
        if gene_name:
            
            gene_id = gene_name.group(1)
        if ref_gene_id:
            
            ref_gene_id = ref_gene_id.group(1)

            mapping[ gene_id ] = '\t'+ref_gene_id+'\t'+"Ref"
            for each_id in all_ids:
                total_mapping[ each_id  ] = ref_gene_id
        else:
            for each_id in all_ids:
                total_mapping[ each_id  ] = raw_id
            mapping[ gene_id ] = "\t\tNew"
END = open( sys.argv[2] , 'w'  )
            
for key,val in mapping.items():
    END.write( key+'\t'+val+'\n' )
    
for key,val in total_mapping.items():
    MAPPING_ALL.write( key+'\t'+val+'\n' )