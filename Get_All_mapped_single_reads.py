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
                       
                          default='Read2',)
	
	
	
	parser.add_option("-m", "--Mapping1", 
	                  action="store",
	                  dest="Mapping1",

	                  help="Read1 Mapping result"
	                  )			





	parser.add_option("-M", "--Mapping2", 
                          dest="Mapping2",
                          help="Read2 Mapping result")

	parser.add_option("-o", "--Output", 
                          dest="Output",
                          help="Output")
	

	(options, args) = parser.parse_args()

	READ1 = fastq_check( open( options.Read1 ,'rU')  )
	READ2 = fastq_check( open( options.Read2 ,'rU')  )
	MAPPING1 = fastq_check( open( options.Mapping1,'rU'  ) )
	MAPPING2 = fastq_check( open( options.Mapping2,'rU'  ) )
	OUTPUT1 = open(options.Output+'.read1.fastq','w')
	OUTPUT2 = open(options.Output+'.read2.fastq','w')
	
	'''得到所有的没有mapping上的reads'''
	all_need_hash = {}
	for each_data in [MAPPING1, MAPPING2  ]:
		for a,b,c,d in each_data:
			title = a[1:-2]
			all_need_hash[ title ] = ''
	for a,b,c,d in READ1:
		title = a[1:-2]
		if title in all_need_hash:
			OUTPUT1.write(a+b+c+d)
	for a,b,c,d in READ2:
			title = a[1:-2]
			if title in all_need_hash:
				OUTPUT2.write(a+b+c+d)	