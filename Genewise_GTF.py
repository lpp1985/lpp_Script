#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/9/5
"""

from lpp import *
from copy import deepcopy



if __name__ == '__main__':
	RAW = open(sys.argv[1])
	data_l = re.split("\n{2,}", RAW.read())
	END = open(sys.argv[2], 'w')
	for e_b in data_l:
		if not e_b.strip():
			continue
		exon_list = []
		cds_list = []
		for line in e_b.split("\n"):
			line_l = line.split("\t")
			if line_l[2] == "mRNA":
				line_l[-2] = "0"
				end_line = ""
				trans_id = re.search("ID\=([^\;]+)", line).group(1) + '_t'
				gene_id = re.search("Parent\=([^\;]+)", line).group(1) + '_g'
				new_l =  deepcopy(line_l)
				if line_l[6] == '+':
					new_l[2] = "start_codon"
					new_l[4] = str(int(line_l[3]) + 2)
					END.write("\t".join(new_l[:-1]) + '\tgene_id "%s"; transcript_id "%s";\n' % (gene_id, trans_id))
					
					new_l[2] = "stop_codon"
					new_l[3] = str(int(line_l[4])  -2)
					new_l[4] = line_l[4]				
					end_line = "\t".join(new_l[:-1]) + '\tgene_id "%s"; transcript_id "%s";\n' % (gene_id, trans_id)
					
				if line_l[6] == "-":
					
					new_l[2] = "stop_codon"
					new_l[4] = str(int(line_l[3]) + 2)
					END.write("\t".join(new_l[:-1]) + '\tgene_id "%s"; transcript_id "%s";\n' % (gene_id, trans_id))
				
					new_l[2] = "start_codon"
					new_l[3] = str(int(line_l[4])-2)
					new_l[4] = line_l[4]

					
					end_line = "\t".join(new_l[:-1]) + '\tgene_id "%s"; transcript_id "%s";\n' % (gene_id, trans_id)
			if line_l[2] in ["CDS", "exon"]:
				line_l[-1] = '\tgene_id "%s"; transcript_id "%s";\n'% (gene_id, trans_id)
				if line_l[2] == "CDS":
					cds_list.append("\t".join( line_l))
				else:
					line_l[-2] = "."
					exon_list.append("\t".join( line_l))
		exon_list = iter(exon_list)
		cds_list = iter(cds_list)
		for key in cds_list:
			END.write(key)
			END.write(exon_list.next())
		END.write(end_line + '\n\n' )