#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/6/1
from lpp import *
import glob
def N50(fasta):
	pool = []
	RAW = fasta_check( open( fasta,'rU' ) )
	distribution=Ddict()
	total = 0
	for t,s in RAW:
		s= re.sub( '\s+','',s )
		length = len( s )
		pool.append(length)
		total +=length
		distribution[ length ][ t ] = ''
	
	def med(newlist):
		cache = 0
		for key in newlist:
			if cache < 0.5*all_sum:
				cache+=key
			else:
				return key

	def N90(  newlist  ):
		cache = 0
		for key in newlist:
			if cache < 0.9*all_sum:
				cache+=key
			else:
				return key
	newlist = []
	pool =sorted(pool)[::-1]
	#for i in pool:
	#	assert i>0 ; assert isinstance( i,int )
	#	newlist+=[i]*i
	newlist = pool
	all_sum = sum(  newlist )
	out_median = med( newlist )
	out_n90 = N90( newlist  )
	out_max = max( distribution )
	out_min = min( distribution )
	base_name = fasta.split('/')[-1].rsplit(".",1)[0]
	END = open(base_name+'.stats','w')
	END.write( 'N50\tN90\tmax\tmin\tTotal\n' )
	END.write( '%s\t%s\t%s\t%s\t%s\n'%( out_median,out_n90, out_max, out_min,total ) )
	DIS = open( base_name+'.dis','w' )
	DIS.write( 'Base\tNO\n' )
	for i in sorted( distribution ):
		DIS.write( '%s\t%s\n'%( i, len(distribution[i]) ) )
	return out_median,distribution,out_max,out_min
fil_l = sys.argv[1:]
for fil in fil_l:
	N50( fil )

