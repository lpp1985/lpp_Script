#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/22
from lpp import *
from optparse import OptionParser
usage = '''usage: python2.7 %prog -i [CLUSTER] -o [Prefix] -f [ FASTA  ]'''
parser = OptionParser(usage =usage ) 
parser.add_option("-i", "--INPUT", action="store", 
                  dest="input", 
                  help="CLUSTER")
parser.add_option("-o", "--prefix", action="store", 
                  dest="prefix", 
                  default = 'Cluster',
                  help="The PREFIX you want!!!!")
parser.add_option("-f", "--fasta", action="store", 
                  dest="fasta", 
                  default = 'Fasta',
                  help="The Fasta file !!!!")
(options, args) = parser.parse_args() 
intput = options.input
PREFIX = options.prefix
FASTA =  fasta_check(  open(   options.fasta,'rU'  )  )
title_seq = {}
for t ,  s in FASTA:
	title_seq[  t[ 1:-1 ].split()[0]  ] =  s

if not os.path.isdir(  PREFIX ):
	os.mkdir(  PREFIX  )
CLUSTER = fasta_check( open(  options.input  )   )
for t,s in CLUSTER:
	END = open (  PREFIX+'/CLUSTER'+re.search( 'Cluster (\d+)',t ).group(1)+'.fasta','w' )
	all_seq = re.findall( '>([^\.]+)' ,s )
	for each_t in all_seq:
		END.write( '>%s\n%s'%(  each_t, title_seq[ each_t ]    ) )
