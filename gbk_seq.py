#!/usr/bin/env python
#coding:utf-8
# Author:   --<>
# Purpose: 
# Created: 2013/1/8

from Bio import SeqIO
from lpp import *
from Bio.SeqRecord import SeqRecord
import itertools

RECORD = SeqIO.parse( sys.argv[1],'genbank' )
GENE = open( sys.argv[2]+'.fa','w' )
from Bio import Seq
SeqRecord.format
i=0
for each_data in RECORD:
    
    seq = each_data.seq
    
    seq = re.sub( '(\w{70})','\\1\n',str(seq)  )
    GENE.write(">CTG%s\n"%(i))	
    GENE.write(seq+'\n')
    i+=1
