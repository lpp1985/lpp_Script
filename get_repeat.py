#!/usr/bin/env python
#coding:utf-8
# Author:   --<>
# Purpose: 
# Created: 2012/9/1

from lpp import *
INPUT = fasta_check(  open( sys.argv[1] ,'rU' )  )
seq = re.sub( '\s+','', INPUT.next(  )[-1] )
LOCATION = open(  sys.argv[2],'rU'  )
END = open( 'repeat.fasta','w' )
i=0
for line in LOCATION:
   i+=1 
   print( line )
   start,end = [ int(x.strip()) for x in line.split('\t')     ]
   end_seq = seq[ start:end+1   ]
   END.write( '>%s\n%s\n'%(  i,re.sub(  '(\w{60})','\\1\n',end_seq   )       )  )
   