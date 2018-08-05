#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/7/13
"""

from lpp import *

if __name__=="__main__":


	usage = '''usage: python2.7 %prog'''
	parser = OptionParser(usage =usage ) 

	parser.add_option("-f", "--Fasta", action="store", 
	                  dest="Fasta", 

	                  help="All Protein Fasta Sequence")
	
	parser.add_option("-i", "--Indel", action="store", 
	                  dest="Indel", 
	                  help="Indel Input")
	
	parser.add_option("-s", "--Snp", action="store", 
	                  dest="Snp", 
	                  help="Snp Input")
	
	parser.add_option("-k", "--SnpResult", action="store", 
	                  dest="SnpResult", 
	                  help="Snp Result")
	

	parser.add_option("-m", "--IndelResult", action="store", 
	                  dest="IndelResult", 
	                  help="Indel Result")
	
	parser.add_option("-l", "--FrameShiftResult", action="store", 
	                  dest="FrameShiftResult", 
	                  help="FrameShift Result")	
(options, args) = parser.parse_args() 



ALL_FASTA = fasta_check( open(options.Fasta))
all_seq = {}
for t, s in ALL_FASTA:
	name = t.split()[0][1:]
	all_seq[name] = ""
	
SNP = open(options.Snp, 'rU')
INDEL = open(options.Indel, 'rU')
SNP_END = open(options.SnpResult, 'w')
INDEL_END = open(options.IndelResult, 'w')
FRAMESHIFT_END =  open(options.FrameShiftResult, 'w')
FRAMESHIFT_END.write("Gene\tseq id\tPosition\tmutation\tannotation\tdescription\n")
SNP_END.write("Gene\tseq id\tPosition\tmutation\tannotation\tdescription\n")

INDEL_END.write( "Gene\tseq id\tMissing start\tMiss End\tMissingSize\n")


gene_ProChange = Ddict()
gene_frameshift = Ddict()
gene_missing = Ddict()
for line in SNP:
	if "evidence" in line or "\t" not in line:
		
		continue
	line_l = line.split("\t")
	if "intergenic" in line_l[4]:
		continue
	gen_name = re.sub("[^\w]+", "", line_l[5])

	result_data = "\t".join([line_l[0], line_l[1], line_l[2], line_l[3], line_l[4], line_l[-1] ])
	if "coding" in line_l[4]:
		gene_frameshift[gen_name] [result_data] = ""
		
	else:
		gene_ProChange[ gen_name] [ result_data] = ""
		
for line in  INDEL:
	if re.search("^\s+", line) or "Unassigned missing coverage evidence" in  line:
		
		continue
	
	line_l = line.split("\t")
	line_res = line_l[3:-4]
	if "intergenic" in line_l[-2]:
		
		continue
	data_all = line_l[-2] + ', ' + line_l[-1][:-1]
	data_all = data_all.replace("[", "").replace("]", "").replace('/', ", ")
	all_gene = set( data_all.replace(", ", " ").split())

	for each_gene in all_gene:
		
		if each_gene in all_seq:
			gene_missing[each_gene]["\t".join(line_res)] = ""
	
for each_gene in gene_missing:
	
	if each_gene in gene_ProChange:
		del gene_ProChange[each_gene]

for each_gene in gene_frameshift:
	
	if each_gene in gene_ProChange:
		del gene_ProChange[each_gene]
	
for each_gene in gene_ProChange:

	for each_snp in  gene_ProChange[each_gene ]:
		
		SNP_END.write(each_gene + '\t' + each_snp)
	
for each_gene in gene_frameshift:
	for each_snp in  gene_frameshift[each_gene]:
		
		FRAMESHIFT_END.write(each_gene + '\t' +each_snp )

for each_gene in gene_missing:
	for each_snp in  gene_missing[each_gene]:
		
		INDEL_END.write(each_gene + '\t' + each_snp + '\n')
