#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/3/16
"""

import glob
from optparse import OptionParser
usage = "python2.7 %prog [options]"
parser = OptionParser(usage =usage )
parser.add_option("-i", "--input_path", action="store",
                  dest="input_path",

                  help="QC stats folder")

parser.add_option("-o", "--out", action="store",
                  dest="output",

                  help="output")





if __name__ == '__main__':
	(options, args) = parser.parse_args()
	all_data = glob.glob( options.input_path+'/*.total.stats')
	END = open( options.output,'w' )
	RAW = open(all_data[0],'rU')
	END.write(RAW.next())
	RAW.close()
	for f in all_data:
		RAW = open(f)
		RAW.next()
		END.write(RAW.next())
		RAW.close()
	
	