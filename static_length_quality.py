#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/4/29
from lpp import *
from optparse import OptionParser
usage = '''usage: python2.7 %prog -f [Fasta_file] -q [Quality_file] -o [ OUTPUT  ]'''
parser = OptionParser(usage =usage ) 
parser.add_option("-f", "--FASTA", action="store", 
	dest="fasta", 
	help="FASTA FILE")

parser.add_option("-o", "--PREFIX", action="store", 
	dest="prefix", 

	help="The PREFIX you want!!!!")

parser.add_option("-q", "--QUALITY", action="store", 
	dest="quality", 

	help="The quality file!!!!")
parser.add_option("-l", "--LENGTH", action="store", 
	dest="length", 

	help="THE LENGTH  threshold")
(options, args) = parser.parse_args() 


FASTA = fasta_check(open( options.fasta,'rU'  ) )
QUAL = fasta_check(open( options.quality,'rU'  ) )
lenghth = int(  options.length   )
length_hash = Ddict()
qual_hash = []
qual_100_hash = []
average_quality = Ddict()
reads_length   = 0
for t,s in FASTA:
	reads_length+=1
	s = re.sub( '\s+','',s  )
	q = [int(x) for x in QUAL.next()[-1].split()]
	ave_q = sum( q )/len(q)
	average_quality[ ave_q ][ s ] = ''
	length_hash[ len(s) ] [  t ] = ''
	qual_hash+=q
	if len(s)>lenghth:
		qual_100_hash+=q
prefix = options.prefix
TOTAL_LENGTH = open( prefix+'.toal_length','w' )
for s in sorted(  length_hash  ):
	TOTAL_LENGTH.write( '%s\t%s\n'%( s, len( length_hash[s]  )   )   )

AVE = open(  prefix+'.ave_qual','w' )
AVE.write( 'Qual\tNumber\n' )
a30 = len( [ key2 for x in    average_quality for key2 in average_quality[ x ]   if x >=30 ] )
a25 = len( [ key2 for x in    average_quality  for key2 in average_quality[ x ]  if x >=25 and x<30 ] )
other = reads_length - a30 - a25
AVE.write( 'Q30\t%s\n'%( a30 ) )
AVE.write( 'Q25\t%s\n'%( a25 ) )
AVE.write( 'Other\t%s\n'%( other ) )
TOTAL_QUAL = 	open( prefix+'.toal_qual','w' )
TOTAL_QUAL.write( 'Qual\tNumber\n' )
def check_qual( qual_hash, TOTAL_QUAL ):
	q_40 = [x for x in qual_hash if x>=40]
	q_30 = [x for x in qual_hash if x>=30 and x <40]
	q_20 = [x for x in qual_hash if x>=20 and x <30]
	TOTAL_QUAL.write( 'Q40\t%s\n'%( len(  q_40  ) )  )
	TOTAL_QUAL.write( 'Q30\t%s\n'%( len(  q_30  ) )  )
	TOTAL_QUAL.write( 'Q20\t%s\n'%( len(  q_20  ) )  )
	TOTAL_QUAL.write( 'Other\t%s\n'%( len( qual_hash  )- len(  q_40  ) - len(q_30 ) - len(  q_20  ) )  )
check_qual( qual_hash, TOTAL_QUAL )
TOTAL_QUAL_100 = 	open( prefix+'.100_toal_qual','w' )
TOTAL_QUAL_100.write( 'Qual\tNumber\n' )
check_qual( qual_100_hash, TOTAL_QUAL_100 )