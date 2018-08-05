#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/9
'''usage quality_stats.py [path]  [appendix]'''
from lpp import *
'''usage: quality_stats.py [reads] [ prefix ]'''
RAW = fasta_check( open( sys.argv[1],'rU'  ) )
seq = re.sub( '\s+','',RAW.next()[-1] )
all_hash = {}
for i in xrange( 0,len( seq   )-25 ):
	each_frag = seq[i:i+25]
	for frag in [ each_frag,complement( each_frag  )     ]:
		if frag in all_hash:
			all_hash[ frag  ] +=1
		else:
			all_hash[ frag  ] =1
end_hash = {}

for each_k in all_hash:
	time = all_hash[ each_k ]
	if time in end_hash:
		end_hash[ time  ]+=1
	else:
		end_hash[ time  ] =1
END = open( 'Kmer.stats','w'  )
END.write( 'times\tnumber\n'   )
for each_t in sorted(  end_hash ):
	END.write( '%s\t%s\n'%( each_t,end_hash[ each_t ]   )  )