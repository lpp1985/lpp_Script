#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/6/9
from lpp import *
REF_FASTA = fasta_check( open( sys.argv[1] ,'rU')  )
END_FASTA = open( sys.argv[2],'w' )
title = sys.argv[3]
END_FASTA.write(  '>'+title+'\n' )

all_seq = ''
for t,s in REF_FASTA:
	all_seq+=s
all_seq = re.sub('\s+','',all_seq)
all_seq = re.sub( '(\w{60})','\\1\n',all_seq  )
END_FASTA.write( all_seq+'\n' )
