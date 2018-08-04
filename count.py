#!/usr/bin/env python
# coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/25
from lpp import *
# usage :count.py [  bowtie_analysis-end   ][  gene_reads_analysis_count-end  ] [   multiple_reads_end     ]
import sys
contig_reads = Ddict()
reads_contig = Ddict()
for line in open(  sys.argv[1]  ):
	line_l =line.split('\t')
	reads = line_l[1].split(';')
	contig = line_l[0]
	for read in reads:
		contig_reads[  contig ][ read ] = ''
		reads_contig[  read  ][  contig   ] = ''
END = open(sys.argv[2],'w')
OTHER = open(sys.argv[3],'w')
for contig in contig_reads:
	contig_num = 0.0
	for read in contig_reads[ contig ]:
		if len( reads_contig[ read ]  )>1:
			#contig_num+= 1.0/ len(  reads_contig[ read ]  )
			contig_num+= 1.0
			OTHER.write( contig+'\t'+ read+'\t%.4f'%( 1.0/ len(  reads_contig[ read ]  ) )+'\n'  )
		else:
			contig_num+= 1.0
	END .write( contig+'\t'+'%s\n'%( contig_num  )  )