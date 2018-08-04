#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/22
from lpp import *
from optparse import OptionParser
import subprocess
import multiprocessing
usage = '''usage: python2.7 %prog -a [ACE_file] -o [OUTPUT] -t [ Threshold  ]'''
parser = OptionParser(usage =usage ) 

parser.add_option("-a", "--ACE", action="store", 
	dest="ace", 
	help="ace_file")

parser.add_option("-o", "--output", action="store", 
	dest="output", 
	default = 'Output_file',
	help="The PREFIX you want!!!!")

parser.add_option("-t", "--threshold", action="store", 
	dest="threshold", 
    type='int', 
    default = 200, 
	help="The bp number threshold in the fasta  you want !!!!")

parser.add_option("-s", "--singleTon", action="store", 
	dest="singleTon", 


	help="The bp number singletopn in the fasta !!!!")

(options, args) = parser.parse_args() 
ace = options.ace
output = options.output
threshold = options.threshold
singleTon = options.singleTon

os.system(  'ace2contig g -i %s  -o %s  '% (  ace , 'ace.cache'   ) )
number = 0
CACHE = re.split(    '\#\#.+\n',open(  'ace.cache'  ).read()  )

END = open( output,'w')

have = {}

for each_block in CACHE:
	s1 = re.split(  '\#',each_block )[0]
	s1 = re.sub( '\s+','',s1 )
	if len(s1)>=threshold:
		number +=1
		s1 = re.sub( '(\w{60})','\\1\n',s1 )
		END.write( '>%s\n%s\n'%( number,s1 )   )
		
os.remove(  'ace.cache'  )

SINGLETON = fasta_check( open( singleTon,'rU'  ) )
for t,s in SINGLETON:
	s1 = re.sub( '\s+','',s )
	if len(s1) >=threshold:
		number += 1
		END.write( '>%s\n%s\n'%( number,s1 )   )

