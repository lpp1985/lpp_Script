#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose:
# Created: 2011/8/29
from lpp import *
import multiprocessing
from optparse import OptionParser
from quality_stats_singleEnd  import *
'# you could type [  SCRIPT_NAME  ] -h to see the help log !!!!'

# get the sample name from file name
def get_sample( name  ) :
	sample = os.path.split(os.path.abspath(name)  )[-1].split('.')[0]
	return sample



class fastq_quality_pair(  object    ):

	## input a hash of a copy of  locals----input_hash = locals().copy(), the set the files variant to be your input file_list, just like input_hash[ 'files'  ] = [ read1,read2     ]
	def __init__(  self,**kword  ):

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

	def next( self  ):
		def check( output_hash ):
			for key in [ 'Filter','failFilter','Trim','failTrim'    ]:
				if key in output_hash:
					status1 = key
					return status1
		## the status_transfer_hash can make output more accuracy , and output_hash record the output file to be output
		status_transfer_hash = { 'Filter':'Trim' ,'failFilter' :'Trim','Trim':'Trim','failTrim':'failTrim'  }
		File_hash = { 'Filter':[ self.FILTERPAIR1,self.FILTERPAIR2 ], }
		if self.trim :
			File_hash [ 'Trim'] = [ self.TRIMPAIR1, self.TRIMPAIR2      ] 
			
		try:
			output_hash1 = self.read1.next()
			output_hash2 = self.read2.next()
		except StopIteration:
			try:
				output_hash2 = self.read2.next()

			except StopIteration:
				subdevide = self.read1.subdevide
				STATS = open(  self.read1.stats_path+ self.SampleName+'.total.stats' ,'w' )
				STATS.write('#Sample\tAllRawData(MB)\tAllRawRead\tAllFilterData(MB)\tAllFilterRead\tDataFilterPerc(%)\tReadFilterPerc(%)\tN(%)\tQ20(%)\tGC(%)\n')
				self.RawData = self.read1.raw_data_count + self.read2.raw_data_count

				self.RawReads = self.read1.raw_reads_count + self.read2.raw_reads_count

				self.FilterReads =  self.filter_pair_read

				self.FilterData =  self.filter_pair_data/self.read1.subdevide

				self.FilterPerc =  100*self.FilterData/ self.RawData

				self.FilterReadPerc =  100*self.FilterReads/ float(  self.RawReads  )

				self.N_perc = 100* (  self.read1.N+self.read2.N  ) / (    subdevide *float(   self.read1.raw_data_count +self.read2.raw_data_count  )     )   
				
				self.GC_perc = 100*( self.read1.GC+ self.read2.GC  )/ (  subdevide*float( self.read1.raw_data_count + self.read2.raw_data_count )   )
				self.Q20_perc = 100*(self.read1.Q20+self.read2.Q20)/(  subdevide*float(   self.read1.raw_data_count + self.read2.raw_data_count  ) )

				STATS.write( self.SampleName+'\t'+ \
							 '\t'.join([  isinstance( x,int )  and   '%s'%(  x)  or '%.2f'%( x ) for x in   \
								        [self.RawData,
								         self.RawReads,
								         self.FilterData,
								         self.FilterReads ,
								         self.FilterPerc,
								         self.FilterReadPerc ,
				                         self.N_perc,
				                         self.Q20_perc,
				                         self.GC_perc ]  \
								        ]) +'\n'\
							 )
				raise StopIteration







		#self.log.write(  output_hash1[ 'pairname'   ]  )

		#if output_hash1[ 'pairname'   ] != output_hash2[ 'pairname'   ]:
			#message = '%s in read1, %s in read2 of line %s seems not to be a pair of reads,please check!!!!!'%( output_hash1[ 'pairname'   ] ,  output_hash2[ 'pairname'   ],  self.read1.linenumber-3)
			#raise ValueError,message
		status1 = check(  output_hash1  )

		status2 = check(  output_hash2  )

		if status1== status2:
			if status1 =='Filter':
				self.filter_pair_data += len(  output_hash1[  status1 ].split('\n')[1]  )
				self.filter_pair_data += len(  output_hash2[  status2 ].split('\n')[1]  )
				self.filter_pair_read +=2

			if status1 in File_hash:
				File_hash[ status1  ][0].write( output_hash1[  status1 ]   )
				File_hash[ status1  ][-1].write( output_hash2[  status2 ]   )
		else:
			status_transfer1 = status_transfer_hash[ status1 ]

			status_transfer2 = status_transfer_hash[ status2 ]

			if status_transfer1==status_transfer2:

				if status_transfer1 in File_hash:
					File_hash[ status_transfer1  ][0].write( output_hash1[  status_transfer1 ]   )
					File_hash[ status_transfer1  ][-1].write( output_hash2[  status_transfer2 ]   )
		#self.log.write( '\t'+status1+'\t'+status2+'\n'   )

def PairEndRun(  ):

	def quality_check( input_files    ):

		files = all_file_group[ input_files ]
		data = rootparamater
		data.update( globals() )
		data.update(  locals()  )

		aa = fastq_quality_pair( **data )
		for i in aa :
			pass

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
	#Input is ->files[  two path tuple  ],  trim[ boolen ], quality [ int], threshold[ float >0 ]  , trim_length[ int, >0  ] 
	QCPAIR_class = fastq_quality_pair( **paramater  )
	for x in QCPAIR_class:
		pass
	#the output variance is --->  QCPAIR_class  , raw_reads( int,  1 ) , raw_data( float, Mb ), filter_reads(  int,1 ) , filter_data(  float,Mb  ), filter_perc(  float, 100*perc ), N_perc(  float, 100*perc ), GC_perc(  float, 100*perc ),Q20_perc(  float,100*perc ) 
	return QCPAIR_class,QCPAIR_class.RawReads , QCPAIR_class.RawData ,  QCPAIR_class.FilterReads  , QCPAIR_class.FilterData, QCPAIR_class.FilterPerc, QCPAIR_class.N_perc , QCPAIR_class.GC_perc, QCPAIR_class.Q20_perc


if __name__ == '__main__':
	path,all_raw_file,quality,threshold,trim,trim_length = get_paramater()
	PairEndRun()
