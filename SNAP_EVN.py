#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/8/16
"""

from lpp import * 



RAW = fasta_check( open(sys.argv[1], 'rU'))
END = open(sys.argv[2], 'w')
for t, s in RAW:
	
	seq = t[1:-1]
	name = ""
	for line in s.split("\n"):
		line_l = line.split("\t")
		if line_l[0] == "Esngl":
			name = line_l[-1].strip()
			END.write( seq + '\tSNAP\tgene\t' + line_l[1] + '\t' + line_l[2] + '\t.\t' + line_l[3] + '\t.' + '\t' + "ID=%s;Name=%s\n" % (name, name) )
			END.write( seq + '\tSNAP\tmRNA\t' + line_l[1] + '\t' + line_l[2] + '\t.\t' + line_l[3] + '\t.' + '\t' + "ID=%s_mRNA;Parent=%s\n" % (name, name) )
			END.write( seq + '\tSNAP\texon\t' + line_l[1] + '\t' + line_l[2] + '\t.\t' + line_l[3] + '\t.' +  '\t' + "ID=%s_exon;Parent=%s_mRNA\n" % (name, name) )
			END.write( seq + '\tSNAP\tCDS\t' + line_l[1] + '\t' + line_l[2] + '\t.\t' + line_l[3] + '\t' + line_l[7] + '\t' + "ID=%s_cds;Parent=%s_mRNA\n" % (name, name) )
		elif   line_l[0] == "Einit":
			i = 1
			cache = ""
			start = line_l[1]


			cache += seq + '\tSNAP\texon\t' + line_l[1] + '\t' + line_l[2] + '\t.\t' + line_l[3] + '\t.' +  '\t' + "ID=%s_exon%s;Parent=%s_mRNA\n" % (name, i,name) 
			cache += seq + '\tSNAP\tCDS\t' + line_l[1] + '\t' + line_l[2] + '\t.\t' + line_l[3] + '\t' + line_l[7] + '\t' + "ID=%s_cds%s;Parent=%s_mRNA\n" % (name, i, name) 
			i+=1
		elif  line_l[0] == "Exon":
			cache += seq + '\tSNAP\texon\t' + line_l[1] + '\t' + line_l[2] + '\t.\t' + line_l[3] + '\t.' +  '\t' + "ID=%s_exon%s;Parent=%s_mRNA\n" % (name, i,name) 
			cache += seq + '\tSNAP\tCDS\t' + line_l[1] + '\t' + line_l[2] + '\t.\t' + line_l[3] + '\t' + line_l[7] + '\t' + "ID=%s_cds%s;Parent=%s_mRNA\n" % (name, i, name) 
			i+=1
		elif line_l[0] == "Eterm":
			
			end = line_l[1]
			cache =seq + '\tSNAP\tmRNA\t' + start + '\t' + end + '\t.\t' + line_l[3] + '\t.' + '\t' + "ID=%s_mRNA;Parent=%s\n" % (name, name) + cache
			cache = seq + '\tSNAP\tgene\t' + start + '\t' + end + '\t.\t' + line_l[3] + '\t.' + '\t' + "ID=%s;Name=%s\n" % (name, name)  + cache
			cache += seq + '\tSNAP\texon\t' + line_l[1] + '\t' + line_l[2] + '\t.\t' + line_l[3] + '\t.' +  '\t' + "ID=%s_exon%s;Parent=%s_mRNA\n" % (name, i,name) 
			cache += seq + '\tSNAP\tCDS\t' + line_l[1] + '\t' + line_l[2] + '\t.\t' + line_l[3] + '\t' + line_l[7] + '\t' + "ID=%s_cds%s;Parent=%s_mRNA\n" % (name, i, name) 
			END.write(cache)
