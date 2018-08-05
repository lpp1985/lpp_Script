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
for line in RAW:
    line_l = line.split("\t")
    line_l[0] = line_l[0].split("|")[0]
    new_data = line_l[:2]    
    if "\ttranscript\t" in line:
        new_data.append("gene")
        new_data.extend([line_l[3],line_l[4],'.',line_l[6],'.'])
        gene_id = re.search("gene_id \"(\S+)\"\;",line).group(1)
        transcript_id=re.search("transcript_id \"(\S+)\"\;",line).group(1)
        attribute = "ID=%s;Name=%s"%( gene_id,gene_id )
        END.write('\t'.join(new_data)+'\t'+attribute+'\n')
        new_data[2] = "mRNA"
        transcript_id=re.search("transcript_id \"(\S+)\"\;",line).group(1)
        attribute = "ID=%s;Parent=%s"%( transcript_id,gene_id )
        END.write('\t'.join(new_data)+'\t'+attribute+'\n')
    else:
        new_data.append("exon")
        new_data.extend([line_l[3],line_l[4],'.',line_l[6],'.'])
        exon_number= re.search("exon_number \"(\d+)\"",line).group(1)
        attribute = "ID=%s.exon%s;Parent=%s"%( transcript_id,exon_number,transcript_id )
        END.write('\t'.join(new_data)+'\t'+attribute+'\n')
        new_data[2] = "CDS"
        attribute = "ID=cds.%s;Parent=%s"%( transcript_id,transcript_id )
        END.write('\t'.join(new_data)+'\t'+attribute+'\n')
        
        