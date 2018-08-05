#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/5/9
"""
import re,os
from optparse import OptionParser

	
if __name__ == "__main__":
	usage = '''usage: python2.7 %prog [options] 
'''
	parser = OptionParser(usage =usage )

	parser.add_option("-i", "--Input", action="store",
	                  dest="inputData",
	                  help="Input Data")	
	parser.add_option("-o", "--Output", action="store",
	                  dest="Output",
	                  help=" Output")		
	parser.add_option("-c", "--Caption", action="store",
	                  dest="Caption",
	                  help=" Caption")			
	(options, args) = parser.parse_args()
	inputData = os.path.abspath( options.inputData )
	Output = os.path.abspath( options.Output )
	Caption = options.Caption
	OUTPUT = open(  Output,'w' )
	
	OUTPUT.write("""
\\begin{figure}[H]

    \\centering
    \\includegraphics[width=0.8\\textwidth]{%s}
    \\captionsetup{labelsep=period}
    \\caption{%s}
\\end{figure}
	
	
	"""%(inputData,Caption))
		
	
