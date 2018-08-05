#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose:QC singleEND on Each Node
# Company: Chinese Human Genomic Center at Beijing
# Created: 2011/12/30
from lpp import *
from itertools import imap
from optparse import OptionParser
import copy_reg
import os
import multiprocessing
from types import MethodType

# The two kind of error to report
class CouldNotWrite(  IOError ):
	pass




class fastq_quality_class(  fastq_check    ):
	
	phred_64_quality_hash={}
	for x in xrange(0,64):
		phred_64_quality_hash[ chr( x+33  ) ] = x
	@staticmethod
	def check_path( total_path , check= False ):


		if check:
			'''Check the abspath, if it does not exist , It will automatically report error'''
			[ path, fil ]=os.path.split(  os.path.abspath( total_path ) )
			if os.path.exists( path ) and os.path.isfile( total_path ):
				return [ path, fil ]
			else:
				error_message = 'Path %s does not exist!!'%(path)
				raise FileNotExistError, error_message
		try:
			os.makedirs( total_path  )
		except:
			pass


		return os.path.abspath( total_path  )+os.sep
	def __init__( self,**kword  ):
		self.subdevide =1048576
		self.Queue = multiprocessing.Queue()
		self.RAWHANDLE  = open( os.path.abspath( kword['input_data'] ), 'rU' )

		self.RAW = fastq_check( self.RAWHANDLE    )
		self.filter_read_count =0
		self.filter_data_count = 0.0
		self.raw_reads_count=0
		self.raw_data_count = 0.0
		[ root_path, file_name ]  = self.check_path( kword[ 'input_data'  ] ,check = True    )
		self.block_size = kword[ 'block_size'  ]
		self.cache_size = kword[  'cache_size'  ]
		self.file_name = file_name
		self.sample_name  = file_name.split('.')[0]



		self.filter_path = root_path+'/Filter/%s/'%( self.sample_name )
		self.check_path(  self.filter_path )
		self.GC = 0
		self.trim = kword['trim']
		self.N = 0
		self.Q20 = 0
		self.stats_path = root_path+'/stats/'
		self.check_path(  self.stats_path  )		

		## the quality hash to represent the quality in each site,which can make searching more fast

		## self.attribute record the variant input,quality ->quality threshold, threshold -> to threshold you
		##want, length -> reads_length after trim to recieve, fail_file -> low quality reads, trim-> True: to
		##trim, False: Not trim, succ_trim -> trim successfuly file to record, fail_trim -> trim failed reads to
		##record, output-> file which filtered reads write in
		self.quality = kword['quality']
		isinstance( self.quality,int )
		self.threshold = kword['threshold']
		isinstance( self.threshold ,(float,int) )

		# the filter fail reads will write to this file
		try:
			self.FAIL_FILE = open( self.filter_path+ file_name+ '.failFilter' ,'w')
		except IOError:
			message = '%s could not write, Please check the permission of that path!!'%(  self.filter_path+ file_name+ '.failFilter'  )
			raise CouldNotWrite, message
		if  self.trim :

			self.trim_path = root_path+'/Trim/%s/'%( self.sample_name )


			self.check_path(  self.trim_path   )
			self.trim_data_count = 0.0
			self.trim_read_count = 0
			self.trim = True
			self.length = kword[ 'trim_length'  ]
			try:
				self.SUCCTRIM = open( self.trim_path + file_name+'.Trim'  ,'w')
				self.FAILTRIM = open( self.trim_path + file_name +'.failTrim' ,'w')
			except IOError:
				message = '%s could not write, Please check the permission of that path!!'%(  self.trim_path  )
				raise CouldNotWrite, message
		try:
			self.OUTPUT = open( self.filter_path+  file_name +'.Filter', 'w'    )
		except IOError:
			message = '%s could not write, Please check the permission of that path!!'%(  self.filter_path )
			raise CouldNotWrite ,message
		## decide the threshold is locate in (0,1), then it should represent low quality percentage
		## else it represent the total number of low quality base insequence,the thre_status is 'number'
		## represents this
		## situation else its value is none

	## This function's input is  quality line of fastq and output is True of False of this reads could be reserved
	## quality is the quality threshold to filter, quality_line in the quaility_line of fastq, low_quality is
	## the low_quality base number of fastq,low_threshold is the threshold variant to record threshold , if
	## thre-status == number, it represents filter threshold is low quality base number else it is low quality
	## percentage
	@staticmethod
	def quality_plot( QUE, FILE_HANDLE    ):
		try:
			box_graph = '%s.qualstats.png'%( FILE_HANDLE.name )
			distr_graph = '%s.nucdistr.png'%(  FILE_HANDLE.name )
			os.system(  '''( fastx_quality_stats -i %s > %s.qualstats &&'''%(  FILE_HANDLE.name,FILE_HANDLE.name  ) + '''fastx_quality_boxplot.R %s.qualstats %s  2>&1  >/dev/null &&'''%(  FILE_HANDLE.name,box_graph   ) + '''fastx_nucleotide_distributionPer.R  %s.qualstats  %s  2>&1  >/dev/null )'''%(  FILE_HANDLE.name,distr_graph   ) )
			QUE.put(  [box_graph , distr_graph] )
		except:
			raise AttributeError,'The %s is plot Fail'%(  FILE_HANDLE.name )
	
	
	def mutiprocess_plot(  self  ):
		RAWQUEUE = multiprocessing.Queue(  )
		FILTERQUEUE = multiprocessing.Queue(  )
		RAWPLOT_PROCESS = multiprocessing.Process(target=self.quality_plot ,args = ( RAWQUEUE , self.RAWHANDLE, ) )
		FILTERPLOT_PROCESS = multiprocessing.Process(target=self.quality_plot ,args = ( FILTERQUEUE , self.OUTPUT ,) )
		FILTERPLOT_PROCESS.start()
		RAWPLOT_PROCESS.start()
		FILTERPLOT_PROCESS.join()
		RAWPLOT_PROCESS.join()
		
		[self.FilterBoxGraph, self.FilterDisGraph] =FILTERQUEUE.get()
		[self.RawBoxGraph, self.RawDisGraph] = RAWQUEUE.get()
	
		
	
	
	def threshold_check( self,quality_line   ):

		quality_line = re.sub( '\s+','' ,quality_line)
		total_length = len( quality_line )
		low_quality = len(  [ x     for x in quality_line if self.phred_64_quality_hash[x]<self.quality    ] ) ## this is the base number of  total low quality base

		if self.threshold > 1:
			low_threshold = low_quality
		else:
			low_perc = float( low_quality  ) / total_length
			low_threshold = low_perc
		if low_threshold > self.threshold:
			return False
		else:

			return True

	@classmethod
	def static( cls,  each_block  ):

		GC  = len(    re.findall( '([G|C])',each_block[1]    )     )

		N   = len(  re.findall( '(N)',each_block[1]     )    )

		Q20 = len(  [ x     for x in re.sub( '\s+$','',  each_block[-1] ) if cls.phred_64_quality_hash[x]>=20    ] )

		Data= len( each_block[1][:-1] )
		return GC, N, Q20, Data

	#a function to do qc
	def combiner(self  ):
		Queue = self.Queue
		
		filteredData = ''	
		trimedData =''
		failfilteredData = ''
		failtrimedData = ''				
		while Queue.qsize():	
				
			data = self.Queue.get()
			self.raw_reads_count+= data[0]
			self.raw_data_count +=data[1]
			if self.trim:
				self.trim_read_count += data[2]
				self.trim_data_count += data[3]
			self.filter_read_count += data[4]
			self.filter_data_count += data[5]

			self.Q20 += data[6]
			self.GC += data[7]
			self.N += data[8]
			filteredData += data[9]
			trimedData += data[10]
			failfilteredData+= data[11]
			failtrimedData += data[12]
	
		self.OUTPUT.write( filteredData  )
		
		self.FAIL_FILE.write(  failfilteredData )
		if self.trim:
			self.SUCCTRIM.write(  trimedData )
			self.FAILTRIM.write(  failtrimedData )		
	def trimQC(  self,each_block ):

		name,seq,define2,quality = each_block


		seq_re = re.search( '^N*([^N]\S{%s,}[^N])N*\n$'%( self.length-2 ) ,seq     )

		if  seq_re:
			[ trim_start,trim_stop ]     = seq_re.span(1)

			seq_trim     = seq[ trim_start:trim_stop  ]
			quality_trim = quality[ trim_start : trim_stop  ]
			output_block = [ name , seq_trim+'\n' , define2,quality_trim+'\n' ]



			out_hash['Trim_GC']  , out_hash['Trim_N']  ,out_hash['Trim_Q20'] ,out_hash['TrimData']  = self.static( output_block    )



			status = 'Trim'

			out_hash[ status ] = output_block




		else:
			status = 'failTrim'

			out_hash[ status ] =  each_block



		return (   out_hash     )




	def filterQC( self,each_block  ):


		quality  = each_block[-1]
		seq = each_block[1]
		if self.threshold_check(  quality   ):

			status = 'Filter'
			out_hash[ status ] = each_block
			out_hash['Filter_GC']  , out_hash['Filter_N']  ,out_hash['Filter_Q20'] ,out_hash['FilterData']  = self.static( each_block    )



		else:
			status = 'failFilter'
			out_hash[ status ] = each_block
		return out_hash





	def Multiple_QC( self,each_block  ):
		global out_hash
		out_hash = {}

		[name_line,seq,define2,quality ] = each_block
		out_hash['RAW_GC']  , out_hash['RAW_N']  ,out_hash['RAW_Q20'] ,out_hash['RAW_Data']  = self.static( each_block    )
		if self.trim:
			self.trimQC( each_block  ) 
			if 'Trim' in out_hash:
				self.filterQC(  out_hash[  'Trim' ]  ) 
		else:

			self.filterQC(  each_block  ) 
		return out_hash

	
	
	def each_qc(  self, each_iter  ):
						
		def getall( list_data ):
			return ''.join(list_data )			
		def check( each_hash ):
			global GC,N,Q20,raw_data_count,raw_reads_count_iner,filter_read_count,filter_data_count,trim_read_count,trim_data_count,filter_detail_cache
			GC +=each_hash[ 'RAW_GC' ] 
			Q20 += each_hash[ 'RAW_Q20' ] 
			
			N +=each_hash['RAW_N'] 

			raw_data_count += each_hash['RAW_Data']
			raw_reads_count_iner += 1					
			
			
			
			if  'Filter' in each_hash:
				filter_cache[ each_hash[ 'Filter' ] [0 ] [:-3] ] =  each_hash[ 'Filter' ]  
				filter_read_count+=1
				filter_data_count += each_hash[ 'FilterData' ] 
				filter_detail_cache[ each_hash[ 'Filter' ] [0 ] [:-3] ]  = each_hash
			elif 'failFilter'  in each_hash :

				fail_filter_cache.append( ''.join(each_hash[ 'failFilter' ] ) )	
			
			if  'Trim' in each_hash:
				trim_cache[   each_hash[  'Trim'  ] [ 0 ] [:-3]   ] = each_hash[ 'Trim' ]
				trim_read_count +=1
				trim_data_count += each_hash[ 'TrimData' ] 	
			elif 'failTrim' in each_hash:
				
				fail_trim_cache.append(  ''.join(each_hash[ 'failTrim' ]  )  )
						
							
													
				
		
		
		
					
					
		filter_cache ={}
		trim_cache ={}
		fail_filter_cache = []
		fail_trim_cache = []
		
		global GC,N,Q20,raw_data_count,raw_reads_count_iner,filter_read_count,filter_data_count,trim_read_count,trim_data_count,filter_detail_cache
		filter_detail_cache = {}
		GC = 0
		Q20 = 0
		N = 0
		raw_data_count = 0.0
		raw_reads_count_iner = 0			
		filter_read_count = 0
		filter_data_count = 0.0
		trim_read_count = 0
		trim_data_count = 0.0				
		
		map(   check  , each_iter)

		all_filterData, all_trimData =  map(  getall,[     [  ''.join(y)   for x,y in filter_cache.items()  ] ,       [  ''.join(y)   for x,y in trim_cache.items()  ]    ] )
		all_fail_trimData , all_fail_filterData= map( getall,    [fail_trim_cache, fail_filter_cache  ]  )

		self.Queue.put( [  
                                                        raw_reads_count_iner,
                                                        raw_data_count,
                                                        trim_read_count,
                                                        trim_data_count,
                                                        filter_read_count,
                                                        filter_data_count,
                                                        Q20,
                                                        GC,
                                                        N,
                                                        all_filterData,
                                                        all_trimData,
                                                        all_fail_filterData,
                                                        all_fail_trimData,
                                                        
                                                        
                                                                   ]      )	
		
		return   filter_detail_cache	
	def outputStatic(self):
		
		subdevide = self.subdevide

		STATS = open( self.stats_path + re.sub( '(\.[^\.]+)$','', self.file_name   ) +'.stats'  ,'w')
		[name,category ]= self.file_name.split('.')[:2]

		self.N_perc = 100*self.N/float(self.raw_data_count )
		self.GC_perc  = 100*self.GC /float(   self.raw_data_count    )
		self.Q20_perc = 100*self.Q20/float(   self.raw_data_count   )		
		self.raw_data_count_out    = float(  self.raw_data_count  ) /subdevide
		self.filter_data_count_out = float(  self.filter_data_count )/subdevide
		self.filter_read_perc    = 100 * self.filter_read_count/float( self.raw_reads_count )
		self.filter_data_perc    = 100 * self.filter_data_count/float(   self.raw_data_count   )

		



		if self.trim:
			self.trim_data_count_out   = float(  self.trim_data_count )/subdevide
			self.trim_read_perc    = 100 * self.trim_read_count/float( self.raw_reads_count)
			self.trim_data_perc    = 100 * self.trim_data_count/self.raw_data_count
			self.stats_title = '#Sample\tCategory\tRawReads\tRawData(MB)\tTrimReads\tTrimData(MB)\tTrimReadPerc(%)\tTrimDataPerc(%)\tFilterReads\tFilterData(MB)\tFilterReadsPerc(%)\tFilterDataPerc(%)\tN(%)\tGC(%)\tQ20(%)'
			output_rank  = [
	                        self.raw_reads_count,
	                        self.raw_data_count_out  ,
	                        self.trim_read_count,
	                        self.trim_data_count_out,
	                        self.trim_read_perc ,
	                        self.trim_data_perc,
	                        self.filter_read_count,
	                        self.filter_data_count_out,
	                        self.filter_read_perc,
	                        self.filter_data_perc,
	                        self.N_perc,
	                        self.GC_perc,
	                        self.Q20_perc
	                ]
		else:
			self.stats_title = '#Sample\tCategory\tRawReads\tRawData(MB)\tFilterReads\tFilterData(MB)\tFilterReadsPerc(%)\tFilterDataPerc(%)\tN(%)\tGC(%)\tQ20(%)'
			output_rank  = [

	                        self.raw_reads_count,
	                        self.raw_data_count_out,
	                        self.filter_read_count,
	                        self.filter_data_count_out,
	                        self.filter_read_perc,
	                        self.filter_data_perc,
	                        self.N_perc,
	                        self.GC_perc,
	                        self.Q20_perc
	                ]


		output = []
		for each_data in output_rank:
			if isinstance( each_data, int ):

				output.append( '%s'%( each_data )  )
			else:
				output.append( '%.2f'%( each_data )  )

		output_data = '\t'.join( output   )

		STATS.write(  self.stats_title+'\n'  )
		STATS.write( name+'\t'+category+'\t'+output_data+'\n' )			
		self.mutiprocess_plot()		
	def Qc(  self ):
				
							
				
			
		for out_iter in self:

			process = multiprocessing.Process(    target = self.each_qc, args=[ out_iter,  ]             )
			process.start()
			
			while len(multiprocessing.active_children() ) >self.cache_size or  self.Queue.qsize()>self.cache_size:
				
				self.combiner(   )
				

			
			

		else :

			while len(multiprocessing.active_children() ) :
				self.combiner(   )
			self.outputStatic()
			






	def __iter__(self):
		return self
	#give the result of the data
	def next(self):


		#the block number of now,if equal to block number then run qc,else put data init until enough of end
		block_count = 0
		#cache_list is a list to store cache data
		cache_hash = {}
		#out_hash store the output result
		out_hash = {}

		status = ''
		try:
			while block_count<self.block_size:



				cache_hash[   tuple(self.RAW.next())       ] = ''
				block_count+=1

		except StopIteration:
			pass
		if not cache_hash:
			raise StopIteration
		out_iter = imap( self.Multiple_QC ,cache_hash    )
		return out_iter










	def get_WebNeed( self ):
			#the output variance is ---> raw_reads( int,  1 ) , raw_data( float, Mb ), filter_reads(  int,1 ) , filter_data(  float,Mb  ), filter_perc(  float, 100*perc ), N_perc(  float, 100*perc ), GC_perc(  float, 100*perc ),Q20_perc(  float,100*perc ), FilterFilePath( char,  ), RawBoxGraph( fileName ),  RawDisGraph( fileName ) ,FilterBoxGraph( fileName ), FilterDisGraph( fileName )
		
			return self.raw_reads_count , self.raw_data_count, self.filter_read_count ,self.filter_data_count,self.filter_data_perc, self.N_perc, self.GC_perc, self.Q20_perc, self.OUTPUT, self.RawBoxGraph, self.RawDisGraph, self.FilterBoxGraph, self.FilterDisGraph	
# get the paramater of commandline
def get_paramater(  ):
	try:
		'# you could type [  SCRIPT_NAME  ] -h to see the help log !!!!'
		usage='''usage: python %prog [options]

	    It can stats quality and filter bad quality reads'''
		parser = OptionParser(usage =usage )
		parser.add_option("-p", "--Path", action="store",
		                  dest="path",
		                  type='string',
		                  help="Input path")
		parser.add_option("-a", "--Appendix", action="store",
		                  dest="appendix",
		                  type='string',
		                  default='fastq',
		                  help="Input FASTQ")
		
		parser.add_option("-b", "--Block", action="store",
				                  dest="blocks",
				                  type='int',
				                  default=50000,
				                  help="each block size to map reduce")	
		
		
		
		parser.add_option("-c", "--Cache", action="store",
				                                  dest="cache",
				                                  type='int',
				                                  default=10,
				                                  help="The cache size of to Output")			





		parser.add_option("-t", "--Trim", action="store_true", default= False,
		                  dest="trim",
		                  help="Do you want to trim the fastq file( trim the N in head or in tail  )")


		parser.add_option("-l", "--Trim_length", action="store",
		                  dest="trim_length",
		                  type= "int",
		                  default=60,
		                  help="The length of threshold  to  filter after trim, THIS PARAMATER OUGHT TO BE INGORANCE IN [-t] IS MISSING!!!")


		parser.add_option("-q", "--Quality", action="store", type = "int",
		                  dest="quality",
		                  default=20,
		                  help="Minimum quality score to keep.")


		parser.add_option("-r", "--Threshold", action="store", type = "float",
		                  dest="threshold",
		                  default = 0.2,
		                  help='''The Threashold you want to filter, if this paramater's value locate in the area of (0,1), It represent the maxium percentage the low quality accounts for ,
			      if Integer it represents the maxium number of low quality base in the whole reads to trim''')

		(options, args) = parser.parse_args()

		path = options.path
		block_size = options.blocks
		cache_size =  options.cache
		## Split file into different group

		all_raw_file = glob.glob( os.path.abspath(path)+'/*.'+ options.appendix)







		trim = options.trim

		trim_length = options.trim_length

		quality = options.quality
		threshold = options.threshold
		return path , all_raw_file , quality , threshold , trim , trim_length,block_size,cache_size
	except:
		print(  'The paramater you inpout is wrong !! please use -h to see the help!!!'  )
		sys.exit()


def SingleEndRun( ):
	# to run qc to each file
	def quality_check(  each_f ):
		input_data = each_f
		data = {   }
		data .update( globals() )

		data.update( locals() )
		QC_class = fastq_quality_class( **data )
		QC_class.Qc()


	other_paramater = locals()
	for each_f in all_raw_file:

		multiprocessing.Process( target =quality_check, args = ( each_f,  )   ).start()


# run the proecss to qc
def SingleDbRun(  **paramater ):
	# the input variance is [  input_data , trim,  quality, trim_length ,threshold       ]
	QC_class = fastq_quality_class( **paramater  )
	QC_class.Qc()
	#the output variance is ---> raw_reads( int,  1 ) , raw_data( float, Mb ), filter_reads(  int,1 ) , filter_data(  float,Mb  ), filter_perc(  float, 100*perc ), N_perc(  float, 100*perc ), GC_perc(  float, 100*perc ),Q20_perc(  float,100*perc ), FilterFilePath( char,  ), RawBoxGraph( fileName ),  RawDisGraph( fileName ) ,FilterBoxGraph( fileName ), FilterDisGraph( fileName )

	return QC_class.raw_reads_count , QC_class.raw_data_count_out, QC_class.filter_read_count ,QC_class.filter_data_count_out,QC_class.filter_data_perc, QC_class.N_perc, QC_class.GC_perc, QC_class.Q20_perc, QC_class.OUTPUT, QC_class.RawBoxGraph, QC_class.RawDisGraph, QC_class.FilterBoxGraph, QC_class.FilterDisGraph,QC_class.SUCCTRIM.name



if __name__ =='__main__':
	path,all_raw_file,quality,threshold,trim , trim_length ,block_size,cache_size= get_paramater()
	#for input_data in all_raw_file:
		
	#aa = fastq_quality_class(  input_data= all_raw_file[0], quality = quality, trim = trim , trim_length = trim_length, threshold = threshold ,block_size =  block_size,cache_size = cache_size                           )
	SingleEndRun()