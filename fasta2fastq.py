#!/usr/bin/python
#coding:utf-8
# Author:   --<>
# Purpose: 
# Created: 2013/10/18
from lpp import *
RAW = fasta_check( open( sys.argv[1] ))
END = open( sys.argv[2],'w'  )
for t,s in RAW:
    t = t[1:-1]
    s = re.sub('\s+','',s ).upper()
    qua = 'I'*len( s )
    END.write( '@'+t+'\n'+s+'\n'+'+'+t+'\n'+qua+'\n'  )
