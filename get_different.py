#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/6/20
'''get_different.py [MRNA]  [PEP] [RAW]'''
from lpp import *
all_diff = File_Ddict(open( sys.argv[1],'rU' )).read(1,2)
MRNA = fasta_check(  open( sys.argv[2],'rU' ))
PEP = fasta_check(  open( sys.argv[3],'rU' ))
RAW = fasta_check(  open( sys.argv[4],'rU' ) )
MRNA_diff =  open( sys.argv[2]+'_diff','w' )
for t,s in MRNA:
	mrna_id = re.search( '>(\d+)',t ).group(1)
	if mrna_id in all_diff:
		MRNA_diff.write(t.split()[0]+'\n'+s)
PEP_diff = open( sys.argv[3]+'_diff','w' )
for t,s in PEP:
	pep_id = re.search( '>(\d+)',t ).group(1)
	if pep_id in all_diff:
		PEP_diff.write(t.split()[0]+'\n'+s)
RAW_diff = open( sys.argv[4]+'_diff','w' )
for t,s in RAW:
	raw_id = re.search( '>(\d+)',t ).group(1)
	if raw_id in all_diff:
		RAW_diff.write( t.split()[0]+'\n'+s )