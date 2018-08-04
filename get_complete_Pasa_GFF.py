#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2016/6/28
"""

from lpp import *
RAW = open(sys.argv[1])
all_complete = {}

ALL_CDS= fasta_check(open(sys.argv[2],'rU'))
for t,s in ALL_CDS:
    name = t.split()[0].split("|")[0][1:]
    all_complete[ name ] = ""
END = open(sys.argv[3],'w')
for line in RAW:
    asid = re.findall( "Parent\=([^;\|]+)",line)
    if not asid:
        asid = re.findall( "ID\=([^;\|]+)",line)

    asid = asid[0]

    if asid in all_complete:
        END.write(line)