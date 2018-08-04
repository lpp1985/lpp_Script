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
all_complex = {}
for line in ALL_LOCATION:
	
	start,stop = sorted( [ int(x)  for x in  line.split() ]      )
	COMPLEX.write(  '>Complex%s\n%s\n'%( i ,  s[start:stop  ]   )    )
	
	for x in xrange( start,stop  ):
		all_complex[x] = 'N'
	i+=1
i=0
s2 = ''
for each_seq in s:
	if i in all_complex:
		s2+='N'
	else:
		s2+=each_seq
	i+=1
s_list = re.split( 'N+',s2 )
i=0
END = open('Uncomplex.fasta','w')
for each_seq in s_list:
	END.write( '>Seq%s\n%s\n'%( i,each_seq  ) )
	i+=1
print( 'Complex0+Seq0+Complex1+Seq1+Complex2' )