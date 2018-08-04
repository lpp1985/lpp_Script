#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/4/29
from lpp import *
from optparse import OptionParser
usage = '''usage: python2.7 %prog -f [Fasta_file] -q [Quality_file] -o [ OUTPUT  ]'''
parser = OptionParser(usage =usage ) 
parser.add_option("-f", "--FASTQ", action="store", 
	dest="fasta", 
	help="FASTA FILE")

parser.add_option("-o", "--PREFIX", action="store", 
	dest="prefix", 

	help="The PREFIX you want!!!!")



(options, args) = parser.parse_args() 


FASTQ = fastq_check(open( options.fasta,'rU'  ) )


length_hash = Ddict()
qual_hash = []
average_quality = Ddict()
reads_length   = 0
for title,seq,title2,quality in FASTQ:
	
	t = title[1:-1]
	s = re.sub( '\s+','',seq  )
	reads_length+=len(s)
	q = [ord(x)-64 for x in quality[:-1]]
	for each_q in q:
		if each_q not in average_quality:
			average_quality[each_q]=1
		else:
			average_quality[each_q]+=1

prefix = options.prefix


AVE = open(  prefix+'.qual','w' )
AVE.write( 'Qual\tNumber\tPerc\n' )
a30 = sum( [ average_quality[ x ]  for x in  average_quality   if x >=30 ] )
a20 = sum( [ average_quality[ x ]  for x in  average_quality   if x >=20 ] )
under30_20 =  a20 -a30
other = reads_length - a20
AVE.write( 'Q30\t%s\t%.2f\n'%( a30,float(a30)/reads_length ) )
AVE.write( '20<x<30\t%s\t%.2f\n'%( under30_20, float(under30_20)/reads_length ) )
AVE.write( 'Other\t%s\t%.2f\n'%( other, float(other)/reads_length ) )
AVE.write( 'Q20\t%s\t%.2f\n'%( a20, float( a20  )/reads_length ) )




#TOTAL_QUAL = 	open( prefix+'.toal_qual','w' )
#TOTAL_QUAL.write( 'Qual\tNumber\n' )
#def check_qual( qual_hash, TOTAL_QUAL ):
	#q_40 = [x for x in qual_hash if x>=40]
	#q_30 = [x for x in qual_hash if x>=30 and x <40]
	#q_20 = [x for x in qual_hash if x>=20 and x <30]
	#TOTAL_QUAL.write( 'Q40\t%s\n'%( len(  q_40  ) )  )
	#TOTAL_QUAL.write( 'Q30\t%s\n'%( len(  q_30  ) )  )
	#TOTAL_QUAL.write( 'Q20\t%s\n'%( len(  q_20  ) )  )
	#TOTAL_QUAL.write( 'Other\t%s\n'%( len( qual_hash  )- len(  q_40  ) - len(q_30 ) - len(  q_20  ) )  )
#check_qual( qual_hash, TOTAL_QUAL )
#TOTAL_QUAL_100 = 	open( prefix+'.100_toal_qual','w' )
#TOTAL_QUAL_100.write( 'Qual\tNumber\n' )
#check_qual( qual_100_hash, TOTAL_QUAL_100 )