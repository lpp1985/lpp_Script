#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/6/21
from lpp import *
RAW = fasta_check( open(sys.argv[1],'rU' ) )
all_need = {}
for t,s in RAW:
	mr_id = re.search( '^>(\d+)',t ).group(1)
	all_need[  mr_id  ] = ''
RAW = fasta_check( open(sys.argv[2],'rU' ) )
END = open( sys.argv[3],'w' )
for t,s in RAW:
	mr_id = re.search( '^>(\d+)',t ).group(1)
	if mr_id in all_need:
		END.write(t+s)