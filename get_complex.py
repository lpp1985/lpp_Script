#!/usr/bin/env python
#coding:utf-8
# Author:  LPP
# Purpose: 
# Created: 2011/11/7
# Company: Chinese Human Genomic Center at Beijing
from lpp import *

if __name__=='__main__':
	RAW_FILE = fasta_check( open( sys.argv[1],'rU'  )   )
	for t,s in RAW_FILE:
		s = re.sub( '\s+','',s )
		
ALL_LOCATION = open( sys.argv[2],'rU'  )
COMPLEX  = open( 'Complex.fasta' ,'w' )
i=0
for line in ALL_LOCATION:
	i+=1
	start,stop = sorted( [ int(x)  for x in  line.split() ]      )
	COMPLEX.write(  '>%s\n%s\n'%( i ,  s[start:stop  ]   )    )