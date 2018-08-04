#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/4/11
from lpp import *
'''usage extract_pair.py <read1.corr> <read2.corr> output-prefix <read1> <read2>'''
RAW_1 = fasta_check(  open( sys.argv[1],'rU' )  )
RAW_2 = fasta_check(  open( sys.argv[2],'rU' )  )
ALL_read1 = {}
for t,s in RAW_1:
	title = t.split('#')[0][1:]
	if 'no error found' in t:
		ALL_read1[ title ] = ''
Pair_all1 = {}
Pair_all2 = {}
ONLY_2 = {}
ONLY_1 = {}
for t,s in RAW_2:
	if 'no error found' not in t:
		continue
	title = t.split('#')[0][1:]
	if title  in ALL_read1:
		Pair_all2[ title ] = ''
		Pair_all1[ title ] = ALL_read1[ title ]
	else:
		ONLY_2[  title  ] = s
for each_title in ALL_read1:
	if each_title not in Pair_all1:
		ONLY_1[ title  ] = ''
RAW_fq_1 = fastq_check( open( sys.argv[4],'rU' )  )

END1_single = open( sys.argv[3]+'.sing1','w'  )

PAIR_1 = open( sys.argv[3]+'.pair1.fastq','w'  )
for t_seq,seq,t_qual,qual in RAW_fq_1:
	title = t_seq.split('#')[0][1:]
	if title in Pair_all1:
		PAIR_1.write( ''.join(  [t_seq,seq,t_qual,qual]  ) )
	else:
		END1_single.write( ''.join(  [t_seq,seq,t_qual,qual]  ) )
		
RAW_fq_2 = fastq_check( open( sys.argv[5],'rU' )  )

END2_single = open( sys.argv[3]+'.sing2','w'  )

PAIR_2 = open( sys.argv[3]+'.pair2.fastq','w'  )


for t_seq,seq,t_qual,qual in RAW_fq_2:
	title = t_seq.split('#')[0][1:]
	if title in Pair_all2:
		PAIR_2.write( ''.join(  [t_seq,seq,t_qual,qual]  ) )
	else:
		END2_single.write( ''.join(  [t_seq,seq,t_qual,qual]  ) )