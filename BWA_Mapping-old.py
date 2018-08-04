#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/7/27
from lpp import *
import multiprocessing
from optparse import OptionParser
usage = '''usage: python %prog [options] 

It can automaticly run BWA!!'''
parser = OptionParser(usage =usage ) 

parser.add_option("-r", "--ref", action="store", 
                  dest="ref",
                  type='string',
                  help="the reference seq")

(options, args) = parser.parse_args() 
ref = options.ref

# build index

os.system( 'bwa index -a is %s'%(  ref  )  )

def BWA_MAPPING( (path,all_file)  ):
	name = ''
	for each_f in all_file:
		
		if each_f.endswith( '.read1.fastq' ):
			name = each_f.split('.')[0]
			output_preifx = path+'/'+name 
			read1_file = path+'/'+each_f
			#print( ' bwa aln -o 1 -t 30 -q 60 %s %s -f %s_1.sai 1>%s_1.log 2>&1& ' %( ref , read1_file,output_preifx ,output_preifx  ) )
			os.system(  ' bwa aln -o 1 -t 30 -q 60 %s %s -f %s_1.sai 1>%s_1.log 2>&1 ' %( ref , read1_file,output_preifx ,output_preifx  )  )
	for each_f in all_file:
		if each_f.endswith( '.read2.fastq' ):
			name = each_f.split('.')[0]
			read2_file = path+'/'+each_f
			#print(   ' bwa aln -o 1 -t 30 -q 60 %s %s -f %s_2.sai 1>%s_2.log  2>&1& ' %( ref , read2_file  ,  output_preifx , output_preifx )   )
			os.system(  ' bwa aln -o 1 -t 30 -q 60 %s %s -f %s_2.sai 1>%s_2.log  2>&1 ' %( ref , read2_file  ,  output_preifx , output_preifx )   )
		
	if name:
		os.system(  'bwa sampe  -P -N 5 %s %s  %s  %s  %s -f %s.sam 1>%s_PE.log 2>&1'%(  ref , path+'/'+name+'_1.sai'  , path+'/'+name+'_2.sai' ,read1_file ,read2_file ,output_preifx ,output_preifx )    )
		#print(  'bwa sampe  -P -N 5 %s %s  %s  %s  %s -f %s.sam 1>%s_PE.log 2>&1'%(  ref , path+'/'+name+'_1.sai'  , path+'/'+name+'_2.sai' ,read1_file ,read2_file ,output_preifx ,output_preifx )   )


output_list = []
for a,b,c in os.walk(  './'  ):
	output_list.append( [a,c]   )
pool = multiprocessing.Pool( processes=2 ) 
pool.map(  BWA_MAPPING ,   output_list )