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
		if line.startswith("---"):
			break
	trnatuple = namedtuple("Name","scaffold,Number,Begin,End,type,codon,IntronBegin,IntronEnd,Score")._make
	all_has = {}
	for line in RAW:
		data_list =trnatuple([x.strip() for x in line[:-1].split("\t") ])
		begin = data_list.Begin
		end = data_list.End
		intrbeg = data_list.IntronBegin
		intrend = data_list.IntronEnd
		if int(begin)>int(end):
			frame="-"
			begin,end = end,begin
			intrbeg,intrend = intrend,intrbeg
		else:
			frame='+'		
		if data_list.type !="Pseudo":
			if data_list.scaffold in all_has:
				all_has[data_list.scaffold]+=1
			else:
				all_has[data_list.scaffold]=1
			gene_id = data_list.scaffold +'.tRNA.TU.%s'%( all_has[ data_list.scaffold ] )
			mrna_id = data_list.scaffold +'.tRNA.%s'%( all_has[ data_list.scaffold ] )
			GFF.write(
			    "%s\ttRNAscan-SE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n"%(
			        data_list.scaffold,
			        begin,
			        end,
			        frame,
			        gene_id,
			        gene_id
			    )
			           )
			
			GFF.write(
			    "%s\ttRNAscan-SE\ttRNA\t%s\t%s\t.\t%s\t.\tID=%s;Parent=%s;Name=%s;product=%s\n"%(
			        data_list.scaffold,
			        begin,
			        end,
			        frame,
			        mrna_id,
			        gene_id,
			        mrna_id,
			        'tRNA-'+data_list.type
			    )
			)			
			SEQ.write( '>'+mrna_id+'\n'   )
			exon_num = 0
			if int(data_list.IntronBegin) ==0:
				exon_num+=1
				exon_id = mrna_id+'.exon%s'%(exon_num)
				GFF.write(
				    "%s\ttRNAscan-SE\texon\t%s\t%s\t.\t%s\t.\tID=%s;Parent=%s;Name=%s;product=%s\n"%(
				        data_list.scaffold,
				        begin,
				        end,
				        frame,
				        exon_id,
				        mrna_id,
				        exon_id,
				        'tRNA-'+data_list.type
				    )
				)	
				sequence = all_seq[ data_list.scaffold ][ int(begin):int(end)   ]
				SEQ.write(sequence+'\n')
				
			else:
				need_scaffold = all_seq[ data_list.scaffold ]
				sequence = need_scaffold[int(begin):int(intrbeg)]+need_scaffold[int(intrend):int(end)]
				SEQ.write(sequence+'\n')
				exon_num+=1
				exon_id = mrna_id+'.exon%s'%(exon_num)
				GFF.write(
					"%s\ttRNAscan-SE\texon\t%s\t%s\t.\t%s\t.\tID=%s;Parent=%s;Name=%s;product=%s\n"%(
						data_list.scaffold,
						begin,
						intrbeg,
						frame,
						exon_id,
						mrna_id,
						exon_id,
						'tRNA-'+data_list.type
					)
				)					
					
				exon_num+=1
				exon_id = mrna_id+'.exon%s'%(exon_num)
				GFF.write(
					"%s\ttRNAscan-SE\texon\t%s\t%s\t.\t%s\t.\tID=%s;Parent=%s;Name=%s;product=%s\n"%(
						data_list.scaffold,
						intrend,
						end,
						frame,
						exon_id,
						mrna_id,
						exon_id,
						'tRNA-'+data_list.type
					)
				)					
		