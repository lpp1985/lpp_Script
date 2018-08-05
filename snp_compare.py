#!/usr/bin/env python
#coding:utf-8
# Author:   --<>
# Purpose: 
# Created: 2013/6/8

#compare 2nd and 3rd sequencing project and ouput result
#Usage: $Script  2nd.snps 3rd.snps Result
from lpp import *
sorted_genes = {}
def get_data( RAW ):
	data_hash = Ddict()
	RAW.next()
	global sorted_genes
	for line in RAW:
		line_l = line.split('\t')
		if line_l[0]:
			sorted_genes[line_l[0]]=int( line_l[1] )
			name = line_l[0]
		data_hash[name][  int(line_l[1]) ] = '\t'.join( line_l[1:] )
	return data_hash
SEC = open( sys.argv[1],'rU' )
THIRD = open( sys.argv[2],'rU' )
sec_hash = get_data( SEC )
thr_hash = get_data( THIRD )
END = open( sys.argv[3],'w'  )
END.write( 'GENE_ID\tStatus\tGenome Pos\tRef\tAlt\tGene\tSeq\tChange\tAA Change	Category\tLevel\n' )
for each_gene in sorted( sorted_genes,key = lambda x: int( sorted_genes[x]  ) ):
	END.write( each_gene )
	if each_gene in sec_hash and each_gene in thr_hash:
		aa = sec_hash[ each_gene ].keys()
		bb = thr_hash[ each_gene ].keys()  
		aa.extend(bb)
		all_location = set( aa )
		for each_location in sorted( all_location ):
			if each_location in sec_hash[ each_gene ] and each_location in thr_hash[ each_gene ]:
				if sec_hash[ each_gene ][  each_location ] == thr_hash[ each_gene ][  each_location ] :
					status = 'Identical'
					END.write( '\t'+status+'\t'+sec_hash[ each_gene ][  each_location ]   )
				else:
					status = 'Different_2nd'
					END.write( '\t'+status+'\t'+sec_hash[ each_gene ][  each_location ]   )
					status = 'Different_3rd'
					END.write( '\t'+status+'\t'+thr_hash[ each_gene ][  each_location ]   )
			elif each_location in sec_hash[ each_gene ] :
				status = '2nd_unique'
				END.write( '\t'+status+'\t'+sec_hash[ each_gene ][  each_location ]   )
			else:
				status = '3rd_unique'
				END.write( '\t'+status+'\t'+thr_hash[ each_gene ][  each_location ]   )				
					
	elif each_gene in sec_hash:
		status = '2nd_unique'
		for each_location in sorted(  sec_hash[ each_gene ] ):
			END.write( '\t'+status+'\t'+sec_hash[ each_gene ][  each_location ]   )		
	else:
		status = '3rd_unique'
		for each_location in sorted(  thr_hash[ each_gene ] ):
			END.write( '\t'+status+'\t'+thr_hash[ each_gene ][  each_location ]   )		
		
	
