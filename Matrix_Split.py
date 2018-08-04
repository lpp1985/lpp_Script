#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/4/29
"""

from lpp import *
from optparse import OptionParser 
usage = ''' To do Desq of Reads mapping matrix !
Usage python %prog     [Options]''' 
parser = OptionParser(usage =usage ) 
parser.add_option("-o", "--OUTPUTPATH", action="store", 
                  dest="outputpath",
                  default = 'Static_output', 
                  help="OUTPUTPATH")

parser.add_option("-i", "", action="store", 
                  dest="Matrix", 
                  default = 'matrix',
                  help="Matrix of Reads number")

(options, args) = parser.parse_args()

Matrix = options.Matrix

outputpath = options.outputpath+'/'
if not os.path.exists( outputpath  ):
	os.makedirs(outputpath)
RAW = open(Matrix,'rU')
title = RAW.next()
title_l = title.rstrip().split("\t")
gene_express = Ddict()
for line in RAW:
	line_l = line.rstrip().split("\t")
	for i in xrange(1,len(line_l)):
		gene_express[  line_l[0] ][title_l[i]] = line_l[i]

for x in xrange(1,len(title_l[1:])):
	for k in xrange(x+1,len(title_l)):
		END= open(outputpath+title_l[x]+"__"+title_l[k]+'.express','w')
		END.write("gene\t"+'\t'.join( [  title_l[x],title_l[k] ]  )+'\n')
		for key in gene_express:
			END.write(key+'\t'+gene_express[key][ title_l[x] ]+'\t'+gene_express[key][ title_l[k] ]  +'\n'   )
			
