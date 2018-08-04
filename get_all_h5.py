#!/usr/bin/env python
#coding:utf-8
# Author:   --<>
# Purpose: 
# Created: 2013/2/20


from lpp import *
from optparse import OptionParser
if __name__=='__main__':
	usage = '''usage: python2.7 %prog [options] 
	'''
	parser = OptionParser(usage =usage )
	parser.add_option("-p", "--path", action="store",
		              dest="path",
		              default = './',
		              help="The input path you want!!")
	parser.add_option("-o", "--out", action="store",
		              dest="output",
		              default = 'input.fofn',
		              help="The output path  you want!!")
	(options, args) = parser.parse_args()
	path = options.path
	output = options.output	
	END = open( output,'w' )
	for a,b,c in os.walk( path  ):
		abs_path = os.path.abspath(  a  )
		for each_f in c:
			if each_f.endswith( 'bas.h5'  ):
				END.write(   abs_path+'/'+each_f+'\n'   )
				
				