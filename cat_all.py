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
parser.add_option("-i", "--input", action="store", 
                  dest="inputpath",
                  default = './', 
                  help="OUTPUT")




(options, args) = parser.parse_args() 

output   = options.output

appendix = options.append

inputpath = options.inputpath

END = open(output,'w')

seq_all = Ddict()
for a,b,c in os.walk(os.path.abspath(inputpath)):
	for f in c:
		if f.endswith( '.%s'%( appendix ) ):
			try:
				RAW = open(a+'/'+f ,'rU')
	
				for line in RAW:
					END.write( line )
				
			except:
				print( f )


	