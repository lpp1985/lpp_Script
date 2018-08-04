#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/8/29
from lpp import *
import multiprocessing
from optparse import OptionParser
from quality_stats_singleEnd_forStatic import *
'# you could type [  SCRIPT_NAME  ] -h to see the help log !!!!'
usage='''usage: python %prog [options] 

It can stats quality and filter bad quality reads'''

#(options, args) = parser.parse_args() 
#path = options.path

## Split file into different group

all_raw_file = glob.glob( path+'/*.'+ options.appendix)

all_file_group = {}
for each_f in all_raw_file:
	name = re.search( '^([^\.]+)',os.path.split(os.path.abspath(each_f)  )[-1] ).group(1)
	name_suffix = re.escape( name )
	all_file_group[ name] = [  x   for x in all_raw_file if re.search( '^('+name_suffix+')',  os.path.split(os.path.abspath(x)  )[-1]  )     ]





trim = options.trim
if trim :
	status = trim
	trim_length = options.trim_length

quality = options.quality
threshold = options.threshold

class fastq_quality_pair(  object    ):

	## input a hash of a copy of  locals----input_hash = locals().copy(), the set the files variant to be your input file_list, just like input_hash[ 'files'  ] = [ read1,read2     ]
	def __init__(  self,kword  ):

		self.filter_pair_read = 0
		self.filter_pair_data = 0.0
		all_file = sorted(kword[ 'files' ], key= lambda x : int( re.search( 'read(\d)',x,re.I ).group(1) )      )

		self.SampleName = os.path.split( os.path.abspath( all_file[0] ) )[-1] .split('.')[0]


		read1_kword = kword.copy()
		read2_kword = kword.copy()
		read1_kword[ 'input_data'  ] = all_file[0]

		read2_kword[ 'input_data'  ] = all_file[-1]
		self.read1 = fastq_quality_class(  read1_kword  )
		self.read2 = fastq_quality_class(  read2_kword  )
		# this function change the name of raw data name to be paired-name
		def trim_name( char ):
			char = re.sub( '\.read\d.+?$','',char  )
			return char

		#add new name of paired data file
		#if trim:

			#self.read1_trim_pair_output   = open(   trim_name( self.read1.succ_trim.name )+'.Trim.pair1' ,'w' )
			#self.read2_trim_pair_output   = open(   trim_name( self.read2.succ_trim.name )+'.Trim.pair2' ,'w'  )
		#self.read1_filter_pair_output = open(   trim_name( self.read1.output.name )+'.Filter.pair1' ,'w'  )
		#self.read2_filter_pair_output = open(   trim_name( self.read2.output.name )+'.Filter.pair2' ,'w'  )
		#self.log = open(  trim_name ( self.read1.file_handle.name )+'.list','w'    )


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
		#if trim:
			#File_hash = { 'Filter':[self.read1_filter_pair_output,self.read2_filter_pair_output],'Trim':[ self.read1_trim_pair_output, self.read2_trim_pair_output      ] } 
		#else:
			#File_hash = { 'Filter':[self.read1_filter_pair_output,self.read2_filter_pair_output], } 
		try:
			output_hash1 = self.read1.next()
			output_hash2 = self.read2.next()
		except:
			try:
				output_hash2 = self.read2.next()

			except:

				STATS = open(  self.read1.stats_path+ self.SampleName+'.total.stats' ,'w' )
				STATS.write('#Sample\tAllRawData(MB)\tAllRawRead\tAllFilterData(MB)\tAllFilterRead\tDataFilterPerc(%)\tReadFilterPerc(%)\tN(%)\tQ20(%)\tGC(%)\n')
				all_raw_count =( float(self.read1.raw_data_count) + self.read2.raw_data_count  )/1048576

				all_raw_reads = self.read1.raw_reads_count + self.read2.raw_reads_count 

				all_filter_reads =  self.read1.filter_read_count +self.read2.filter_read_count 

				all_filter_data =  ( float(self.read1.filter_data_count) +self.read2.filter_data_count )/1048576

				filter_data_perc =  100*all_filter_data/ float(  all_raw_count  ) 

				filter_reads_perc =  100*all_filter_reads/ float(  all_raw_reads  ) 
				N_perc = 100* (  self.read1.N+self.read2.N  ) /float(   self.read1.raw_data_count +self.read2.raw_data_count  )   ,
				GC_perc = 100*( self.read1.GC+ self.read2.GC  )/float(   self.read1.AT + self.read1.GC +self.read2.GC+ self.read2.AT   ) 
				Q20_perc = 100*(self.read1.Q20+self.read2.Q20)/float(   self.read1.raw_data_count + self.read2.raw_data_count  ) 

				STATS.write( self.SampleName+'\t'+ \
							 '\t'.join(['%.2f'%( x ) for x in   \
								        [all_raw_count, 
								         all_raw_reads,
								         all_filter_data, 
								         all_filter_reads , 
								         filter_data_perc, 
								         filter_reads_perc ,
				                         N_perc, 
				                         Q20_perc, 
				                         GC_perc ]  \
								        ]) +'\n'\
							 )
				raise StopIteration







		#self.log.write(  output_hash1[ 'name'   ]  )

		status1 = check(  output_hash1  )

		status2 = check(  output_hash2  )

		if status1== status2:
			if status1 =='Filter':
				self.filter_pair_data += len(  output_hash1[  status1 ].split('\n')[1]  )
				self.filter_pair_data += len(  output_hash2[  status2 ].split('\n')[1]  )
				self.filter_pair_read +=2

			#if status1 in File_hash:
				#File_hash[ status1  ][0].write( output_hash1[  status1 ]   )
				#File_hash[ status1  ][-1].write( output_hash2[  status2 ]   )
		else:
			status_transfer1 = status_transfer_hash[ status1 ]

			status_transfer2 = status_transfer_hash[ status2 ]

			#if status_transfer1==status_transfer2:

				#if status_transfer1 in File_hash:
					#File_hash[ status_transfer1  ][0].write( output_hash1[  status_transfer1 ]   )
					#File_hash[ status_transfer1  ][-1].write( output_hash2[  status_transfer2 ]   )
		#self.log.write( '\t'+status1+'\t'+status2+'\n'   )




if __name__ == '__main__':
	
	def quality_check( input_files  ):

		files = all_file_group[ input_files ]
		data = globals().copy()
		data.update(  locals()  )
		aa = fastq_quality_pair( data )
		for i in aa :
			pass
	pool = multiprocessing.Pool( processes=20 ) 	
	pool.map( quality_check, all_file_group   )
	#map( quality_check, all_file_group   )
	
	#for each_f in all_file_group:
		#print( each_f )
		#files = all_file_group[ each_f ]
		#aa = fastq_quality_pair( locals() )
		#for i in aa:
			#pass