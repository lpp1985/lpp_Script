#!/usr/bin/python
#coding:utf-8
# Author:   --<>
# Purpose: 
# Created: 2012/11/13

usage = '''usage: python2.7 %prog  -i INPUT -o OUTPUT -l [Length]'''
from lpp import *
from optparse import OptionParser
parser = OptionParser(usage =usage )

parser.add_option("-i", "--INPUT", action="store",
		              dest="ipt",

		              help="The input fasta format file")



parser.add_option("-o", "--OUTPUT", action="store",
                              dest="output",

                              help="The output Result")
parser.add_option("-l", "--length", action="store", 
                  dest="length",
                  default = 0,
                   type= "int",
                  help="length")


(options, args) = parser.parse_args() 

output   = options.output

length =  options.length


ipt = options.ipt
RAW = fasta_check( open( ipt,'rU'   )    )
END = open( output,'w'  )
for t,s in RAW:
	s = re.sub( '\s+','',s  ).upper()
	s = re.sub("\*","N",s)
	if len(s)>= length:
		s = re.sub( '(\w{70})','\\1\n',s  ) 
		if s[:-1]!="\n":
			END.write( t+s+'\n'   ) 
		else:
			END.write( t+s   )

