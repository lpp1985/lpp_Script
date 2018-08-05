#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/20
from lpp import *
from optparse import OptionParser 
usage = '''usage: python2.7 %prog -o OUTPUT -a [Appendix]'''
parser = OptionParser(usage =usage ) 
parser.add_option("-o", "--OUTPUT", action="store", 
                  dest="output",
                  default = 'output', 
                  help="OUTPUT")
parser.add_option("-a", "--append", action="store", 
                  dest="append", 
                  default = 'contig',
                  help="The type you want")

(options, args) = parser.parse_args() 
output   = options.output
appendix = options.append
END = open(output+'.out','w')
seq_all = {}
for a,b,c in os.walk(os.getcwd()):
	for f in c:
		if f.endswith( '.%s'%( appendix ) ):
			
			RAW = open(a+'/'+f ,'rU')
			
			for line in RAW:
				seq = line.split('\t')[-1].replace('\n','')
				if len(seq)>=17:
					seq_all[ line ] = ''
for line in seq_all:
	END.write( line )

	