#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/6/20
'''get_different.py [MRNA]  [PEP] [RAW]'''
from lpp import *
all_diff = File_Ddict(open( sys.argv[1],'rU' )).read(1,2)

RAW = fasta_check(  open( sys.argv[2],'rU' ) )

RAW_diff = open( sys.argv[2]+'_all','w' )
for t,s in RAW:
	raw_id = re.search( '>(\d+)',t ).group(1)
	if raw_id in all_diff:
		RAW_diff.write( t.split()[0]+'\n'+s )