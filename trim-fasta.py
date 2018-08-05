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
parser.add_option("-i", "--input", action="store", 
                  dest="input", 

                  help="The input File")

(options, args) = parser.parse_args() 
output   = options.output
RAW = open(  options.input ,'rU')
END = open(output,'w')
i=1
for line in RAW:
	line_l = line.split('\t')
	END.write( '>%s_'%(i)+line_l[0]+'\n'+line_l[-1]   )
	i+=1

	