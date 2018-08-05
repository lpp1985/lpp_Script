#!/usr/bin/env python
#coding:utf-8
# Author:  LPP
# Purpose: 
# Created: 2011/11/7
from lpp import *

if __name__=='__main__':
	RAW = fasta_check(  open( sys.argv[1],'rU' ) )
	END = open(  sys.argv[2],'w'  )
	quality = sys.argv[3]
	for t,s in RAW:
		END.write( t )
		s = re.sub('\s+','',s       )
		quality_score = [str(quality)]*len(  s  )
		quality_score = '\t'.join(  quality_score    )+'\n'
		END.write( quality_score )
