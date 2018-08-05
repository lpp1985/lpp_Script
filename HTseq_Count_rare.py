#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/11/1
"""

import collections,sys,HTSeq
from optparse import OptionParser
from lpp import *
usage = '''usage: python2.7 %prog [options] Kmer




    Kmer is a list of K value you want,e.g  [ 1, 2, 3, 4 ]'''
parser = OptionParser(usage =usage )



parser.add_option("-s", "--SAM", action="store",
                  dest="sam",
                  help="Sam File")

parser.add_option("-o", "--OUT", action="store",
                  dest="out",
                  help="output")

parser.add_option("-c", "--ReadsCount", action="store",
                  dest="count",
                  help="ReadsCount")

parser.add_option("-r", "--Ref", action="store",
                  dest="ref",
                  help="reference")
(options, args) = parser.parse_args()
sam = options.sam
out = options.out
print(sam)
count = float(options.count)
fasta = fasta_check(open(options.ref,'rU'))
seq_length = {}
for t,s in fasta:
    seq_length[t.split()[0][1:]] = len(re.sub("\s+", "", s))


counts = collections.Counter( )
if sam.endswith("bam"):
    almnt_file = HTSeq.BAM_Reader( sam)
else:
    almnt_file = HTSeq.SAM_Reader( sam)
for almnt in almnt_file:
    if  almnt.aligned and almnt.proper_pair:

        gene_name = almnt.iv
        gene_name=gene_name.chrom

        counts[ gene_name ] += 1

END = open(out,'w')
for gene_id in counts:
    END.write('%s\t%s\n'%(gene_id, 10**11*float(counts[ gene_id ]  )  /seq_length[gene_id]/count )   )
