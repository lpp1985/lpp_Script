#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: Trans 454 to fastq
# Created: 2011/8/8
from lpp import *
from optparse import OptionParser
usage = '''usage: python %prog [options] 

It can automaticly transfer 454 to fastq!!'''
parser = OptionParser(usage =usage ) 
parser.add_option("-s", "--FASTA", action="store", 
                  dest="fasta",
                  type='string',
                  help="the fasta file name")

parser.add_option("-q", "--QUALITY", action="store", 
                  dest="qual",
                  type='string',  
                  help="the quality file name")
parser.add_option("-o", "--OUTPUT", action="store", 
                  dest="output",
                  type='string',  
                  help="THE output fastq file")

(options, args) = parser.parse_args() 


def trans( seq ):
	seq = seq[:-1]
	quality = ''
	for each_number in re.split( '\s+',seq ):
		quality+=chr( int( each_number )+64)  
	return quality
RAW = fasta_check( open( options.fasta ,'rU') )
QUAL = fasta_check( open( options.qual  ,'rU') )
END = open( options.output ,'w' )
for t,s in RAW:
	title = t[1:].replace(' ','__')
	END.write( '@'+title )
	s = re.sub('\s+','',s)
	END.write(s+'\n')
	t,qual_string = QUAL.next()
	quality = trans( qual_string )
	END.write( '+'+ title )
	END.write(quality+'\n')




