#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/11
from lpp import *
fasta = fasta_check(  open( '../../SCF' )  )
seq_name = {}
for t,s in fasta:
	title = t.split()[0][1:] 
	s = re.sub( '\s+','',s )
	seq_name[ title ] = s
LOCI = open( 'combined.loci','rU'  )
FASTA_END = open( 'transcript.fasta','w'  )
MAPPING = open( 'trans_geneName','w'  )
id_trans = {}
for line in LOCI:
	line_l = line[:-1].split('\t')
	transcript_id = line_l[0]

	[( contig,start,end )] = re.findall( '([^\[]+)\[.+?\](\d+)\-(\d+)',line_l[1]   )
	seq = seq_name[  contig  ][ int( start):int( end )     ]
	FASTA_END.write( '>'+transcript_id+'\n'+seq+'\n' )
	all_gene = [x for x in ','.join(line_l[2:]).split(',') if x and x!='-']
	for each_gene  in all_gene:
		id_trans[ each_gene ]  = transcript_id
	MAPPING.write(  '%s\t%s\n'%( transcript_id,'; '.join(  all_gene  )   )  )

for each_f in glob.glob( '*.tmap' ):
	title  = each_f.split('.')[0]
	EXPRE = open( '%s.express'%( title ) ,'w')
	RAW = open( each_f,'rU' )
	RAW.next()
	for line in RAW:
		line_l = line.split('\t')
		cuf_id = line_l[ 4 ]
		fpkm = line_l[ 6 ]
		EXPRE.write(  id_trans[  cuf_id ] +'\t'+fpkm+'\n'  )
		
		
	
	