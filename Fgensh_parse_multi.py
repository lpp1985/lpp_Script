#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/6/2
from lpp import *
RAW = block_reading(  open( sys.argv[1],'rU' )  ,tag = '//')
MRNA = open( sys.argv[1].split('.')[0]+'.mrna' ,'w' )
PROTEIN = open( sys.argv[1].split('.')[0]+'.pep' ,'w' )
DETAIL = open( sys.argv[1].split('.')[0]+'.detail' ,'w' )
ALL_FASTA = fasta_check(  open( sys.argv[2],'rU' )  )
seq_all = {}
for t,s in ALL_FASTA:
	seq_all[t[1:-1]] = re.sub( '\s+','',s )
	
	
for each_b in RAW:
	print( each_b )
	if ' no reliable predictions' in each_b:
		continue
	DETAIL.write( each_b+'\n//\n' )
	
	gene_name = re.search(  'Seq name\: (\S+) ',each_b ).group(1)
	seq_data = each_b.split( 'Predicted protein(s):\n' )[-1]
	pep_list = re.split( '\n(?=\>)', seq_data  )
	for each_gene in pep_list:
		[(start,end)] = re.findall(  )