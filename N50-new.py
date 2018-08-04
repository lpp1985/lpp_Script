#!/usr/bin/python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/7/26
from lpp import *
import glob
def N50(fasta):
	pool = []
	scope_hash = Ddict()
	scope = xrange(0,10000,500)
	for each_scp in scope:
		scope_hash[ each_scp ] = 0
	scope_hash[ '10000' ] = 0

	RAW = fasta_check( open( fasta,'rU' ) )
	
	distribution=Ddict()
	
	for t,s in RAW:
		s= re.sub( '\s+','',s )
		length = len( s )
		pool+=[length]
		distribution[ length ][ t ] = ''
		if length> scope[-1] :
			scope_hash[ '10000' ]+=1
		else:
			for each_scope in scope:
				if length<= each_scope:
					scope_hash[ each_scope ]+=1
					break
		
			
	def median(newlist):
		median_pos = len( newlist )/2
		if len( newlist )%2==0:

			median = float(   newlist[ median_pos ]  +   newlist[ median_pos -1]    )/2
		else:
			median = newlist[ median_pos ]
		return median
	newlist = []
	pool =sorted(pool)
	for i in pool:
		assert i>0 ; assert isinstance( i,int )
		newlist+=[i]*i
	out_median = median( newlist )
	total_number = 0
	out_max = max( distribution )
	out_min = min( distribution )
	for i in distribution:
		for k in distribution[i]:
			total_number+=i
	
	END = open(fasta.split('/')[-1].split(".")[0]+'.stats','w')
	END.write( 'N50\tmax\tmin\tTotal_base\n' )
	END.write( '%s\t%s\t%s\t%s\n'%( out_median, out_max, out_min,total_number  ) )
	DIS = open( fasta.split('/')[-1].split(".")[0]+'.dis','w' )
	DIS.write( 'Base\tNO\n' )
	SCOPE = open( fasta.split('/')[-1].split(".")[0]+'.scope','w' )
	SCOPE.write(  'Range\tNO\n'  )
	for each_scope in scope:
		SCOPE.write( '%s\t%s\n'%( each_scope, scope_hash[ each_scope ]  )  )
	SCOPE.write( '%s\t%s\n'%( '10000', scope_hash[ '10000' ]  )  )
	
	
	for i in sorted( distribution ):
		DIS.write( '%s\t%s\n'%( i, len(distribution[i]) ) )
	return out_median,distribution,out_max,out_min

	
N50( sys.argv[1] )

