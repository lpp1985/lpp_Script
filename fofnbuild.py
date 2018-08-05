#!/usr/bin/python
#coding:utf-8
# Author:   --<>
# Purpose: 
# Created: 2013/8/21

from optparse import OptionParser
from lpp import *

if __name__=='__main__':
	usage = '''usage: python2.7 %prog [options]
	transfer trim overlap relationship'''
	parser = OptionParser(usage =usage )

	parser.add_option("-i", "--INPUT", action="store",
                      dest="inp", 
                      help="Input path") 


	parser.add_option("-s", "--suffix", action="store",
                      dest="suffix",
                      help="suffix of inputfile ,in fasta ,fa,fastq or bas.h5") 

	parser.add_option("-o", "--out", action="store",
                      dest="output",
                      help="OUTPUT Name")


	
	(options, args) = parser.parse_args()
	inp = options.inp
	output = options.output

	suffix = options.suffix	
	OUTPUT = open( output,'w' )
	for a,b,c in os.walk( os.path.abspath( inp  )  ):
		for each_f in c:
			
			if each_f.endswith( suffix ):
				OUTPUT.write( a+'/'+ each_f+'\n'  )
		
		
		
		
		
		
		
	#for each_data in glob.glob(os.path.abspath( inp  )+'/*.'+suffix  ):
		#OUTPUT.write( each_data+'\n'  )
		
