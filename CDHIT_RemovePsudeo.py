#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/9/26
"""

from lpp import * 
from optparse import OptionParser


if __name__ == '__main__':
	usage = '''usage: python2.7 %prog [options]'''
	parser = OptionParser(usage =usage )
	parser.add_option(
	    "-c", "--CLUSTER",
	    action="store",
	    dest="Cluster",
	    help="CD-HIT Cluster list")
	parser.add_option(
	    "-v", "--Validate",
	    action="store",
	    dest="Validate",
	    help="Genewise Validate Gene GFF"	    
	    
	    
	    
	    
	    
	    
	)
	
	
	parser.add_option(
	    "-g", "--GFF",
	    action="store",
	    dest="gff",
	    help="gff result")

	(options, args) = parser.parse_args()
	
	RAW = open( options.Cluster, 'rU')
	all_data = RAW.read()
	all_data_b = re.split("\>Cluster \d+\n", all_data)
	PSUDO_LIST = open("PSUDO.tsv", 'w')
	PSUDO_RELA = open( "PSUDO_Category.tsv", 'w')
	GENE_LIST = open("GENE.tsv", 'w')	
	psudo_dict = {}
	for e_b in all_data_b[1:]:
		representive = re.search("(\d+)nt\, \>(\S+)\.\.\.\s+\*", e_b)
		rep_length = int( representive.group(1) )
		rep_id =  representive.group(2).rsplit(".", 1)[0]
		GENE_LIST.write(rep_id + '\n' )
		all_data = re.findall("(\d+)nt\, \>(\S+)... at", e_b)
		for length, seq_id in all_data:

			seq_id = seq_id.rsplit(".", 1)[0 ]

			length = float(length)
			if length / rep_length < 0.9:
				psudo_dict[seq_id] = ""
				PSUDO_LIST.write(seq_id + '\n')
				PSUDO_RELA.write(seq_id + '\t' +rep_id + '\n' )
			else:
				GENE_LIST.write(seq_id + '\n' )
				
	VALIDATE = open( options.Validate, 'rU')
	cache_data = ""
	for line in VALIDATE:
		if "\tGene\t" in line:
			cache_data += "\n"
		cache_data += line
			
	cache_data =  cache_data.strip() + '\n'
	ALL_GFF = open( "ALL.gff3", 'w')
	GENE_GFF = open( "GENE.gff3", 'w')
	PSUDO_GFF = open( "PSUDO.gff3", 'w')	
	ALL_GFF.write(cache_data)
	GENE_GFF.write(cache_data)
	valiadte_gene = Ddict()
	cache_data = cache_data.strip()
	for e_gene in cache_data.split("\n\n"):
		line_l = e_gene.split("\t")
		scaff = line_l[0]
		start, end = sorted([int(line_l[3]), int(line_l[4])] )
		valiadte_gene[scaff][(start, end)] = ""
		
		
	
	
	
	
	
	
	
	
	
	
	GFF = open( options.gff, 'rU')
	cache_data = ""
	
	for line in GFF:
		if "\tgene\t" in line:
			cache_data += '\n'
		cache_data += line
		
	cache_data = cache_data.strip() + '\n'
	m = 0
	for e_gene in cache_data .split("\n\n"):
		e_gene = e_gene.replace("\tgene\t", "\tGene\t")
		gene_id = re.search("ID\=(g\d+)", e_gene).group(1)
		location = re.search("\tGene\t(\d+)\t(\d+)", e_gene)
		scaff = re.search("(\S+)\tAUGUSTUS\t", e_gene).group(1)
		start = int(location.group(1))
		end =  int(location.group(2))
		start, end = sorted([start, end])
		gene_tag = "no"
		for r_start, r_end in  valiadte_gene[scaff]:
			if end >= r_start and end <= r_end:
				gene_tag = "yes"
				break
			elif   r_end >= start and r_end <= end:
				gene_tag = "yes"
				break
		if gene_tag != "yes":
			
			ALL_GFF.write('\n\n'+e_gene)
			if gene_id in psudo_dict:
				e_gene = e_gene.replace("\tGene\t", "\tPsudo\t")
				PSUDO_GFF.write('\n\n'+e_gene)
			else:
				m+=1
				GENE_GFF.write('\n\n'+e_gene)
				
				
	print(m)
	

	
