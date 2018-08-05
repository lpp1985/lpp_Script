#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/6/19
"""

from lpp import *
from optparse import OptionParser
usage = '''usage: python2.7 %prog [options] Kmer

Convert mauve result to vcf
'''
parser = OptionParser(usage =usage )



parser.add_option("-i", "--INPUT", action="store",
                      dest="input",
                      help="Mauve Input Data")
parser.add_option("-o", "--OUTPUT", action="store",
                  dest="output",
                  help="OutputData")
parser.add_option("-r", "--REFERENCE", action="store",
                  dest="ref",
                  help="Reference_seq")
(options, args) = parser.parse_args()
FASTA = fasta_check(open(options.ref,'rU'))
for t,s in FASTA:
	name = t[1:].split()[0]
END = open(options.output,'w')
END.write("""##fileformat=VCFv4.1
##fileDate=20130606
##reference=NC_000913
##INFO=DP,1,Integer,"Total Depth of Coverage"
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO\n""")
RAW = open(options.input,'rU')
RAW.next()
for line in RAW:
	line_l = line.strip().split("\t")
	ref = line_l[0][0]
	chanded=line_l[0][1:]
	
	END.write('\t'.join(
	   [
	       name,
	       line_l[1],
	       
	       '.',
	       ref,
	       chanded,
	       "255",
	       "PASS",
	       "DP=151"
	       
	   ] 
	)
	          +'\n'
	          )
	
	