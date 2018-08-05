#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/24
from lpp import *
RAW = fasta_check( open(sys.argv[1],'rU')  )
END = open( sys.argv[2],'w' )
for t,s in RAW:
	s1 = re.sub( '\s+','',s  )
	if len(s1)<=200:
		continue
	else:
		END.write( t+s )