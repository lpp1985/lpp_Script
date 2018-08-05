#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: Get taxon catgory add to NR blast annotation
# Created: 2011/5/22
from lpp import *
from optparse import OptionParser
import subprocess
import multiprocessing
usage = '''usage: python2.7 %prog -a [ACE_file] -o [OUTPUT] -t [ Threshold  ]'''
parser = OptionParser(usage =usage ) 

parser.add_option("-n", "--NR", action="store", 
	dest="nr", 
	help="blast file")

parser.add_option("-o", "--output", action="store", 
	dest="output", 
	default = 'Output_file',
	help="The PREFIX you want!!!!")

parser.add_option("-t", "--Taxon", action="store", 
	dest="taxon", 

	help="The taxon gi relationship")


(options, args) = parser.parse_args() 
# Taxon is Taxon file which contain gi_taxon and category like Animal,this file could be producd by taxon_creep.py
nr = options.nr
OUTPUT= open(  options.output ,'w' )
TAXON = open( options.taxon,'rU' )
#NR is the blast annotation file
NR = open( nr,'rU' )
all_have = {}

OUTPUT.write( NR.next()[:-1]+'\tTaxon\n')
gi_annotaton = Ddict()


for line in NR:
	gi = re.search( 'gi\|(\d+)' ,line).group(1)
	gi_annotaton[gi][line[:-1]] = ''
	
for line in TAXON:
	line_l = line.split('\t')
	if line_l[0] in gi_annotaton:
		for each_out in gi_annotaton[ line_l[0] ]:
			OUTPUT.write( each_out+'\t'+line_l[-1] )
			all_have[ line_l[0] ] = ''
			
for each_gi in gi_annotaton:
	if each_gi not in all_have:
		for each_out in gi_annotaton[  each_gi ]:
			print( each_out )


