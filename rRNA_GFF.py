#!/usr/bin/env python
#coding:utf-8
"""
  Author:   ->
  Purpose: 
  Created: 2016/1/11
"""
from lpp import *
from collections import namedtuple

all_seq = {}

if __name__ == '__main__':
	for t,s in fasta_check( open(sys.argv[1],'rU')  ):
		s1 = re.sub("\s+", "", s)
		all_seq[ t[1:].split()[0].split("|")[0]   ] = s1
		
	RAW = open(sys.argv[2],'rU')
	GFF  = open( sys.argv[3],'w')
	SEQ = open(sys.argv[4],'w')
	for line in RAW:
		
		if line.startswith("# --------"):
			break
		
	all_has = {}
	for line in RAW:
		if line.startswith("# --------"):
			break		
		data_list =line.strip().split("\t")
		scaffold = data_list[0]
		
		if scaffold in all_has:
			all_has[scaffold]+=1
		else:
			all_has[scaffold]=1

		gene_id = scaffold +'.rRNA.TU.%s'%( all_has[ scaffold ] )
		
		mrna_id = scaffold +'.rRNA.%s'%( all_has[ scaffold ] )
		data_list[2]="gene"
		data_prefix = "\t".join(data_list[:-1])
		
		GFF.write(
	        "%s\tID=%s;Name=%s\n"%(
		        data_prefix,
	            gene_id,
	            gene_id
	        )
	               )
		data_list[2]="rRNA"
		data_prefix = "\t".join(data_list[:-1])		
		GFF.write(
	        "%s\tID=%s;Parent=%s;Name=%s;product=%s\n"%(
		        data_prefix,
	            mrna_id,
	            gene_id,
	            mrna_id,
	            data_list[-1]
	        )
	    )			
		SEQ.write( '>'+mrna_id+'\n'   )
		exon_num = 0
		
		exon_num+=1
		exon_id = mrna_id+'.exon%s'%(exon_num)
		data_list[2]="exon"
		data_prefix = "\t".join(data_list[:-1])
		GFF.write(
	        "%s\tID=%s;Parent=%s;Name=%s;product=%s\n"%(
		        data_prefix,
	            exon_id,
	            mrna_id,
	            exon_id,
	            data_list[-1]
	        )
	    )	
		sequence = all_seq[ scaffold ][ int(data_list[3]):int( data_list[4])   ]
		if data_list[6]=='-':
			sequence = complement(sequence)
		SEQ.write(sequence+'\n')
			
				
	