#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/4/19
import multiprocessing

from lpp import *
all_file = glob.glob( '*.corr' )
def check( e_f ):
	RAW = fasta_check(open( e_f,'rU' ) )
	cache_hash = {}
	for t,s in RAW:
		cache_hash[ t[1:-1].split(  )[0]   ] = ''
	END = open( e_f.replace( '.coor','.ok' ),'w'   )
	for s_title, seq , q_title,qual in fastq_check( open( e_f.replace( '.coor','' ) ,'rU')   ):
		if s_title[1:-1] in cache_hash:
			END.write( ''.join( s_title, seq , q_title,qual  ) )
thread = int( sys.argv[1] )
pool = multiprocessing.Pool(thread)
pool.map( check, all_file)
print( all_file )