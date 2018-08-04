#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/6/9
from lpp import *
REF1_FASTA = fastq_check( open( sys.argv[1] ,'rU')  )
REF2_FASTA = fastq_check( open( sys.argv[2] ,'rU')  )
END_FASTA = open( sys.argv[3],'w' )
i=1
for a,b,c,d in REF1_FASTA:
	a1 = '>read%s/1\n'%(i)
	a2 = '>read%s/2\n'%(i)
	seq2 = REF2_FASTA.next()[1]
	END_FASTA.write(a1+b+a2+seq2   )
	i+=1
	
