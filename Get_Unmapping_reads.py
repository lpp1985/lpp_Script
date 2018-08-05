#!/usr/bin/env python
#coding:utf-8
# Author:  LPP
# Purpose: 
# Created: 2011/11/7
# Company: Chinese Human Genomic Center at Beijing
from lpp import *
from optparse import OptionParser
if __name__=='__main__':
	usage='''usage: python %prog [options]
	
	It can stats quality and filter bad quality reads'''
	parser = OptionParser(usage =usage )
	parser.add_option("-1", "--Read1", action="store",
                          dest="Read1",
                          type='string',
                          help="Read1")
	parser.add_option("-2", "--Read2", action="store",
                          dest="Read2",
                       
                          default='Read2',
                          help="Input FASTQ")
	
	parser.add_option("-S", "--Single", 
	                  action="store",
	                  dest="Single",

	                  help="Status of Singleend"
	                                  )	
	
	
	
	parser.add_option("-R", "--Reference", 
	                  action="store",
	                  dest="Reference",

	                  help="reference sequence"
	                  )			





	parser.add_option("-O", "--Output", 
                          dest="Output",
                          help="Output result")

	parser.add_option("-M", "--Mapping", 
                          dest="Mapping",
                          help="Already Mapping")
	

	(options, args) = parser.parse_args()

	READ1 = fastq_check( open( options.Read1 ,'rU')  )
	READ2 = fastq_check( open( options.Read2 ,'rU')  )
	SINGLE = open( options.Single,'rU'  )
	Reference = options.Reference
	OUTPUT = open(options.Output,'w')
	
	
	'''得到所有的没有mapping上的reads'''
	single_status = {}
	no_mapping_single = {}
	for line in SINGLE:
		line_l = line.strip().split('\t')
		if line_l[-1] =='R2_have':
			tag = '1'
		else:
			tag = '3'
		no_mapping_single[ line_l[0]+tag   ] = line_l[0]+tag
		single_status[  line_l[0]  ] = line_l[-1]
		
	NOMATCH = open( 'nomatch.fastq','w' )
	for each_direction in [ READ1, READ2    ]:
		for a,b,c,d in each_direction:
			name = a[:-1]
			if name in no_mapping_single:
				NOMATCH.write( a+b+c+d   )
	'''bowtie 开始 mapping'''
	os.system( 'bowtie-build -f %s REF'%( Reference  )  )
	os.system(  'bowtie -c REF -q -n 3 -a -M 1 -p 30 nomatch.fastq >nomatch.list ' )
	MAPPING_END = open( 'nomatch.list','rU' )
	mapping_hash = {}
	for line in MAPPING_END:
		line_l = line[:-1].split('\t')
		title = line_l[0][:-1]
		mapping_hash[  title ] = '\t'.join( line_l[:4]   )
	ALREADYMAPPING = open( options.Mapping ,'rU' )
	for line in ALREADYMAPPING:
		line_l = line.split('\t')
		title = line_l[0][:-1]
		if title in mapping_hash:
			OUTPUT.write( '\t'.join( line_l[:4] )+'\t'+ mapping_hash[ title ]+'\t'+ single_status[ '@'+title    ]+'\n' )