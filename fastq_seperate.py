#!/usr/bin/env python
#coding:utf-8
# Author:  LPP
# Purpose: 
# Created: 2011/11/7
# Company: Chinese Human Genomic Center at Beijing
from lpp import *
from time import strftime
from optparse import OptionParser
import multiprocessing
'# you could type [  SCRIPT_NAME  ] -h to see the help log !!!!'


def sep_data(  DisPath,  BlockNumber ,Inputlist):
	if DisPath[-1] !='/':
		DisPath+='/'	
	sep_QUEUE = multiprocessing.Queue()
	
	timestamp = strftime('%a_%d_%b_%Y_%H_%M_%S/' )	
		
	
	data_name = os.path.split(  Inputlist[0]  )[-1].split(  '.' )[ 0 ]
	timestamp = strftime('%a_%d_%b_%Y_%H_%M_%S/' )
	OutputPath = timestamp+data_name+'/'		
	for i in xrange( 0,BlockNumber  ):
		CachePath = DisPath+'%s/'%(i)+OutputPath
		if not os.path.exists(  CachePath ):
			os.makedirs( CachePath  )	
	def single_sep_data( DisPath,  BlockNumber ,InputData ):	
		'''
		seperate the data into blocks and return the {block:output_path} infomation of data
		'''
		Da_file_name = os.path.split(  InputData )[-1]
			
		def writefunc(  ):
			#write the block 
			for each_tag in cache_hash:
				output_hash[ each_tag].write( cache_hash[  each_tag ]  )
				cache_hash[ each_tag ] = ''	
		output_hash = {}
		cache_hash = {}
		for i in xrange( 0,BlockNumber  ):
			CachePath = DisPath+'%s/'%(i)+OutputPath

			output_hash[i] = open( CachePath +'%s.'%( i )+ Da_file_name ,'w' )
			cache_hash[ i ] = ''
		i=0	
		
		for a,b,c,d in fastq_check( open( InputData,'rU'  )    ):
			i+=1
			tag = i% BlockNumber 
			cache_hash[  tag  ] += a+b+c+d
			
			if i %( 200  )==0:
				writefunc()
		else:
			writefunc()	
		return_hash = {}
		for eac_tag in output_hash:
			return_hash[ eac_tag ] = output_hash[ eac_tag ].name
		sep_QUEUE.put(   return_hash )
	all_process = []	
	for each_data in Inputlist:
		each_process = multiprocessing.Process(      
		        target=single_sep_data,
		        args=[
		             DisPath,BlockNumber,   each_data,
		        ]
		        
		)
		each_process.start()
		all_process.append(  each_process  )
	for each_process in all_process:
		each_process.join()
	all_out = []
	while sep_QUEUE.qsize():
		all_out.append(  sep_QUEUE.get()   )
	return all_out
if __name__=='__main__':
	usage='''usage: python %prog [options]
	
	    It can stats quality and filter bad quality reads'''
	parser = OptionParser(
		usage =usage 
	)

	parser.add_option(
		"-o", "--Onput",
		action="store",
		dest="DisPath",
		type='string',
		help="Distrib path",
		default="../Dis/",
		
		                  )
	
	parser.add_option(
		"-b", "--Block",
		action="store",
		dest="BlockNumber",
		type='int',
		help="Block Number",
		default=10,
		
		                  )
	
	(options, args) = parser.parse_args()
	InputList = args
	DisPath = options.DisPath

	
	inputList = list( set( InputList  ) )
	for 	InputData in InputList:
		if not os.path.isfile( InputData  ) and not os.path.islink( InputData  ) :
			raise IOError,'%s is not exist ,Please Check'%( InputData  )
	BlockNumber = options.BlockNumber
	
	
	sep_data( DisPath,  BlockNumber ,inputList  ) 