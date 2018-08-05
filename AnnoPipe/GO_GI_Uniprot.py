#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2016/7/4
"""
from lpp import *
UNI_GO = open(  sys.argv[1],'rU'   )
UNIPROT = open(  sys.argv[3],'w'   )
UNI_GI = open(  sys.argv[2],'rU'   )
UNI_GORESULT = open(  sys.argv[1]+'.result','w'   )
UNI_GIRESULT = open(  sys.argv[2]+'.result','w'   )
i=0
unip = {}
for line in UNI_GO:
    line_l = line.strip().split("\t")
    if line_l[0] not in unip:
        i= i+1
        unip[ line_l[0] ]=i
        uniId = i
    else:
        uniId = unip[ line_l[0] ]
    UNI_GORESULT.write( str(uniId)+'\t'+line_l[-1]+'\n'  )
    
for line in UNI_GI:
    line_l = line.strip().split("\t")
    gi_list = line_l[-1].split("; ")
    if line_l[0] not in unip:
        i=i+1
        unip[ line_l[0] ]=i
        uniId = i
    else:
        uniId = unip[ line_l[0] ]
    for key in gi_list:
        UNI_GIRESULT.write( str(uniId)+'\t'+key+'\n'  )   
for key,value in unip.items():
    UNIPROT.write( key+'\t'+str(value) +'\n' )
