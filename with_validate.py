#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/6/3
'''usage  with_validate.py [VALIDATE_FASTA] [MATRIX_RAW]  [ MATRIX_END  ]'''
from lpp import *
all_validate = {}
RAW  = fasta_check( open( sys.argv[1],'rU' )  )
for t,s in RAW:
	name = re.search( '\>(\d+)',t  ).group(1)
	all_validate[ name  ] = ''
ALL_MATRIX = open( sys.argv[2],'rU'  )
END_MATRIX = open( sys.argv[3],'w'  )
END_MATRIX.write( ALL_MATRIX.next()  )
print(  len( all_validate )  )
already = {}
for line in ALL_MATRIX:
	line_l = line.split('\t')
	if line_l[0] in all_validate:
		already[ line_l[0] ] = ''
		END_MATRIX.write( line )
#for key1 in all_validate:
	#if key1 not in already:
		#print(key1)