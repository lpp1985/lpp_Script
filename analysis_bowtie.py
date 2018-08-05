#!/usr/bin/python
#Author=LPP
from lpp import *
import sys,re
#analysis_bowtie.py   [bowtie_end_file]  [  output_file   ]
total_mapping = Ddict()
INPUT = open( sys.argv[1] )
OUTPUT = open( sys.argv[2],'w' )
for line in INPUT:
	line_l = line[:-1].split('\t')
	total_mapping[ line_l[2] ][ line_l[0] ] = ''

for contig in total_mapping:
	OUTPUT.write(  contig +'\t'+str( len( total_mapping[ contig ] ) ) +'\n' )
