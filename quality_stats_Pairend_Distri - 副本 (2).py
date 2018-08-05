#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose:
# Created: 2011/8/29
from lpp import *
import multiprocessing
from optparse import OptionParser
from quality_stats_singleEnd_Distri  import *
'# you could type [  SCRIPT_NAME  ] -h to see the help log !!!!'

# get the sample name from file name
def get_sample( name  ) :
	sample = os.path.split(os.path.abspath(name)  )[-1].split('.')[0]
	return sample



class fastq_quality_pair(  object    ):
	
	## input a hash of a copy of  locals----input_hash = locals().copy(), the set the files variant to be your input file_list, just like input_hash[ 'files'  ] = [ read1,read2     ]
	def __init__(  self,**kword  ):
		self.PairQueue = multiprocessing.Queue()
		self.filter_pair_read = 0
		self.filter_pair_data = 0.0
		if not len(   kword[  'files' ] ) ==2:
			print( kword[  'files' ]  )
			message =   'Data %s is error!!! Must only have two files!!!!'%(  '\t'.join(  kword[  'files' ]  ) )
			raise ValueError,message
		try:
			all_file = sorted(   kword[ 'files' ]  ,  key= lambda x : int( re.search( 'read(\d)',x,re.I ).group(1) )      )
		except:
			message = '[%s ]\'s name must fit *.read(\d+).fastq '%(  '\t'.join(  kword[ 'files' ]  )  )
			raise ValueError ,message
		self.SampleName = get_sample( all_file[0]  )


		read1_kword = kword.copy()
		read2_kword = kword.copy()
		read1_kword[ 'input_data'  ] = all_file[0]

		read2_kword[ 'input_data'  ] = all_file[-1]
		self.read1 = fastq_quality_class(  **read1_kword  )

		self.read2 = fastq_quality_class(  **read2_kword  )
		# this function change the name of raw data name to be paired-name
		def trim_name( char ):
			char = re.sub( '\.read\d.+?$','',char  )
			return char

		#add new name of paired data file
		self.trim = kword['trim']
		if self.trim:

			self.TRIMPAIR1   = open(   trim_name( self.read1.SUCCTRIM.name )+'.Trim.pair1' ,'w' )

			self.TRIMPAIR2   = open(   trim_name( self.read2.SUCCTRIM.name )+'.Trim.pair2' ,'w'  )
			self.read2.TRIMPAIR= self.TRIMPAIR2.name
		self.FILTERPAIR1 = open(   trim_name( self.read1.OUTPUT.name   )+'.Filter.pair1' ,'w'  )
		self.FILTERPAIR2 = open(   trim_name( self.read2.OUTPUT.name   )+'.Filter.pair2' ,'w'  )
		self.read1.FILTERPAIR= self.FILTERPAIR1.name
		self.read2.FILTERPAIR= self.FILTERPAIR2.name
		self.log = open(    trim_name ( self.read1.RAWHANDLE.name   )+'.list','w'    )


	def __iter__(self):

		return self

	def eachqc( self ,read1_block, read2_block):
		global filter_pair_read ,filter_pair_data,filter_read1_cont,filter_read2_cont
		filter_pair_read = 0
		filter_pair_data = 0.0
		filter_read1_cont = ''
		filter_read2_cont = ''
		def data_check( paired_read):
			global filter_pair_read ,filter_pair_data,filter_read1_cont,filter_read2_cont
			for each_data in paired_read:
				filter_pair_data += filter_read1[each_data]['FilterData']
				filter_pair_read += 1
				filter_read1_cont += getall( filter_read1[each_data]['Filter'] )
				filter_pair_data += filter_read2[each_data]['FilterData']
				filter_pair_read += 1
				filter_read2_cont += getall( filter_read2[each_data]['Filter']	 )
















		def getall( list_data ):
			return ''.join(list_data )

		filter_read1 = self.read1.each_qc( read1_block )
		filter_read2 = self.read2.each_qc(  read2_block )
		paired_read = filter( lambda x: x in filter_read1, filter_read2    )
		data_check(  paired_read )

		self.PairQueue.put(        [filter_pair_read ,
		                           filter_pair_data ,
		                           filter_read1_cont,
		                           filter_read2_cont ]
		                           )

	def outputStatic( self  ):
		subdevide = self.read1.subdevide
		STATS = open(  self.read1.stats_path+ self.SampleName+'.total.stats' ,'w' )
		STATS.write('#Sample\tAllRawData(MB)\tAllRawRead\tLow_complexity\tAllFilterData(MB)\tAllFilterRead\tDataFilterPerc(%)\tReadFilterPerc(%)\tN(%)\tQ20(%)\tGC(%)\n')
		self.RawData = (self.read1.raw_data_count + self.read2.raw_data_count)/self.read1.subdevide

		self.RawReads = self.read1.raw_reads_count + self.read2.raw_reads_count
		self.low_complex = self.read1.low_complex_data + self.read2.low_complex_data
		self.FilterReads =  self.filter_pair_read

		self.FilterData =  self.filter_pair_data/self.read1.subdevide

		self.FilterPerc =  100*self.FilterData/ self.RawData

		self.FilterReadPerc =  100*self.FilterReads/ float(  self.RawReads  )
		self.low_complexPerc = 100*self.low_complex/ (self.read1.raw_data_count + self.read2.raw_data_count)
		self.N_perc = 100* (  self.read1.N+self.read2.N  ) / (    float(   self.read1.raw_data_count +self.read2.raw_data_count  )     )

		self.GC_perc = 100*( self.read1.GC+ self.read2.GC  )/ (  float( self.read1.raw_data_count + self.read2.raw_data_count )   )
		self.Q20_perc = 100*(self.read1.Q20+self.read2.Q20)/(  float(   self.read1.raw_data_count + self.read2.raw_data_count  ) )

		STATS.write( self.SampleName+'\t'+ \
					 '\t'.join([  isinstance( x,(int,str ))  and   '%s'%(  x)  or '%.2f'%( x ) for x in   \
							[self.RawData,
							 self.RawReads,
		                                         self.low_complexPerc,
							 self.FilterData,
							 self.FilterReads ,
							 self.FilterPerc,
							 self.FilterReadPerc ,
		                                         self.N_perc,
		                                         self.Q20_perc,
		                                         self.GC_perc ]  \
							]) +'\n'\
					 )
		#raise StopIteration
	def combiner( self  ):

		self.filter_read1cache = ''
		self.filter_read2cache = ''
		#print( self.PairQueue.qsize() )
		while self.PairQueue.qsize():

			data = self.PairQueue.get()

			self.filter_pair_read  += data[0]
			self.filter_pair_data  += data[1]
			self.filter_read1cache += data[2]
			self.filter_read2cache += data[3]
		def writeData(  FILEHANDLE,data ):
			FILEHANDLE.write( data )
		process1 = multiprocessing.Process(  target= writeData,args=[   self.FILTERPAIR1,self.filter_read1cache  ]  )
		process2 = multiprocessing.Process(  target= writeData,args=[   self.FILTERPAIR2,self.filter_read2cache  ]  )
		process1.start()
		process2.start()
		process1.join()
		process2.join()
		#writeData(  self.FILTERPAIR1,self.filter_read1cache  )
		#writeData(  self.FILTERPAIR2,self.filter_read2cache  )
		#self.FILTERPAIR1.write( self.filter_read1cache  )
		#self.FILTERPAIR2.write( self.filter_read2cache  )
	def Qc( self ):
		cache_size = self.read1.cache_size
		for read1_block, read2_block in self:
			process = multiprocessing.Process(    target = self.eachqc, args=[   read1_block, read2_block ,  ]             )
			process.start()

			while len(multiprocessing.active_children() ) >cache_size or  self.PairQueue.qsize()>cache_size:

				for each_data in [   self.read1, self.read2 ,self ]:
					each_data.combiner()

		else:
			print(  self.read1.raw_reads_count  )
			while len(multiprocessing.active_children() ) :
				self.read1.combiner()

				self.read2.combiner()

				self.combiner()
			for each_read in [  self.read1,self.read2  ]:
				each_read.outputStatic()
				#process = multiprocessing.Process(    target = each_read.outputStatic, args=[      ]             )
				#process.start()
			self.outputStatic()
			#self.read1.outputStatic()
			#self.read2.outputStatic()
	def next( self  ):

		## the status_transfer_hash can make output more accuracy , and output_hash record the output file to be output


		return self.read1.next() , self.read2.next()



def PairEndRun(  ):

	def quality_check( input_files    ):

		files = all_file_group[ input_files ]
		data = rootparamater
		data.update( globals() )
		data.update(  locals()  )

		aa = fastq_quality_pair( **data )
		aa.Qc()
		#for i in aa :
			#pass

	all_file_group = {}
	for each_f in all_raw_file:
		name = get_sample( each_f )
		name_suffix = re.escape( name )
		all_file_group[ name] = [  x   for x in all_raw_file if re.search( '^('+name_suffix+')',  os.path.split(os.path.abspath(x)  )[-1]  )     ]

	rootparamater = locals()

	for each_group in all_file_group:


		cache = multiprocessing.Process( target = quality_check, args = ( each_group, )  )
		cache.start()


def PairDbRun(  **paramater ):


	QCPAIR_class = fastq_quality_pair( **paramater  )

	QCPAIR_class.Qc()
	#the output variance is --->  QCPAIR_class  , raw_reads( int,  1 ) , raw_data( float, Mb ), filter_reads(  int,1 ) , filter_data(  float,Mb  ), filter_perc(  float, 100*perc ), N_perc(  float, 100*perc ), GC_perc(  float, 100*perc ),Q20_perc(  float,100*perc )
	return QCPAIR_class,QCPAIR_class.RawReads , QCPAIR_class.RawData ,  QCPAIR_class.FilterReads  , QCPAIR_class.FilterData, QCPAIR_class.FilterPerc, QCPAIR_class.N_perc , QCPAIR_class.GC_perc, QCPAIR_class.Q20_perc


if __name__ == '__main__':
	path,all_raw_file,quality,threshold,trim , trim_length ,block_size,cache_size,filtered, tolerate ,kmer  = get_paramater()
	PairEndRun()
