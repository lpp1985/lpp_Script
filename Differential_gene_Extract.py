#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/3/16
"""

from lpp import *
usage = "python2.7 %prog [options]"
parser = OptionParser(usage =usage )
parser.add_option("-i", "--Input", action="store",
                  dest="input_path",

                  help="input_path")
parser.add_option("-r", "--Rpkm", action="store",
                  dest="rpkm",

                  help="rpkm")

parser.add_option("-o", "--out", action="store",
                  dest="output",

                  help="output")

parser.add_option("-a", "--appendix", action="store",
                  dest="appendix",

                  help="appendix")
if __name__ == '__main__':
	(options, args) = parser.parse_args()
	all_has = {}
	for a,b,c in os.walk(  options.input_path ):
		for f in c:
			if f.endswith(  options.appendix  ):
				RAW = open(a+'/'+f)
				RAW.next()
				for line in RAW:
					line_l = line.split("\t")
					all_has[  line_l[0] ]= ""
	END = open(options.output,'w')
	RAW = open( options.rpkm,'rU' )
	END.write( RAW.next() )
	for line in RAW:
		line_l = line.split("\t")
		if line_l[0] in all_has:
			END.write(line)
	