#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose:
# Created: 2011/5/9

from lpp import *
from optparse import OptionParser
import os
import multiprocessing
# The two kind of error to report
class CouldNotWrite(  IOError ):
	pass




class fastq_quality_class(  fastq_check    ):

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
		self.RAWHANDLE  = open( os.path.abspath( kword['input_data'] ), 'rU' )
		self.RAW = super(  fastq_quality_class,self     ).__init__( self.RAWHANDLE    )

		self.filter_read_count =0
		self.filter_data_count = 0.0
		self.raw_reads_count=0
		self.raw_data_count = 0.0
		[ root_path, file_name ]  = self.check_path( kword[ 'input_data'  ] ,check = True    )
		self.file_name = file_name
		self.sample_name  = file_name.split('.')[0]



		self.stats_path = root_path+'/stats/'
		self.check_path(  self.stats_path  )
		self.filter_path = root_path+'/Filter/%s/'%( self.sample_name )
		self.check_path(  self.filter_path )
		self.GC = 0
		self.trim = kword['trim']
		self.N = 0
		self.Q20 = 0


		## the quality hash to represent the quality in each site,which can make searching more fast
		self.phred_64_quality_hash={}
		for x in xrange(1,44):
			self.phred_64_quality_hash[ chr( x+33  ) ] = x
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
				self.fail_trim = open( self.trim_path + file_name +'.failTrim' ,'w')
			except IOError:
				message = '%s could not write, Please check the permission of that path!!'%(  self.trim_path  )
				raise CouldNotWrite, message
			self.stats_title = '#Sample\tCategory\tRawReads\tRawData(MB)\tTrimReads\tTrimData(MB)\tTrimReadPerc(%)\tTrimDataPerc(%)\tFilterReads\tFilterData(MB)\tFilterReadsPerc(%)\tFilterDataPerc(%)\tN(%)\tGC(%)\tQ20(%)'
		else:

			self.stats_title = self.stats_title = '#Sample\tCategory\tRawReads\tRawData(MB)\tFilterReads\tFilterData(MB)\tFilterReadsPerc(%)\tFilterDataPerc(%)\tN(%)\tGC(%)\tQ20(%)'
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
	def threshold_check( cls,quality_line   ):
		quality_line = re.sub( '\s+','' ,quality_line)
		total_length = len( quality_line )
		low_quality = len(  [ x     for x in quality_line if cls.phred_64_quality_hash[x]<cls.quality    ] ) ## this is the base number of  total low quality base

		if cls.threshold > 1:
			low_threshold = low_quality
		else:
			low_perc = float( low_quality  ) / total_length
			low_threshold = low_perc
		if low_threshold > cls.threshold:
			return False
		else:

			return True
## This function plot the quality picture of each run ##
	@staticmethod
	def quality_plot( QUE, FILE_HANDLE    ):
		try:
			box_graph = '%s.qualstats.png'%( FILE_HANDLE.name )
			distr_graph = '%s.nucdistr.png'%(  FILE_HANDLE.name )
			os.system(  '''( fastx_quality_stats -i %s > %s.qualstats &&'''%(  FILE_HANDLE.name,FILE_HANDLE.name  ) + '''fastx_quality_boxplot.R %s.qualstats %s&&'''%(  FILE_HANDLE.name,box_graph   ) + '''fastx_nucleotide_distributionPer.R  %s.qualstats  %s )'''%(  FILE_HANDLE.name,distr_graph   ) )
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


	def next(self):
		## The Try excpt capture the last iteration of programe and output the quality of each_file
		## output_hash is a variant record reads_status of a read,status is variant to record reads status,the value
		## name is the reads name. The value failTrim is trim fail-status, the value Trim is represents could cold Trim
		## the value fail_Filter is low quality Filter , value Filter is could be successfuly filter Reads

		output_hash = {}

		status = ''


		try:
			name,seq,define2,quality = super( fastq_quality_class, self  ).next(   )

			raw_cont =''.join( [  name,seq,define2,quality  ]   )

			self.GC += len(    re.findall( '([G|C])',seq    )     )

			self.N += len(  re.findall( '(N)',seq    )    )

			self.Q20 += len(  [ x     for x in re.sub( '\s+$','',  quality ) if self.phred_64_quality_hash[x]>=20    ] )

			self.raw_data_count+= len( seq[:-1]  )
			self.raw_reads_count +=1

		except StopIteration:
			subdevide = 1048576
			self.subdevide =subdevide

			STATS = open( self.stats_path + re.sub( '(\.[^\.]+)$','', self.file_name   ) +'.stats'  ,'w')
			[name,category ]= self.file_name.split('.')[:2]

			self.raw_data_count    = float(  self.raw_data_count  ) /subdevide
			self.filter_data_count = float(  self.filter_data_count )/subdevide
			self.filter_read_perc    = 100 * self.filter_read_count/float( self.raw_reads_count )
			self.filter_data_perc    = 100 * self.filter_data_count/float(   self.raw_data_count   )
			self.N_perc = 100*self.N/(self.raw_data_count * subdevide )
			self.GC_perc  = 100*self.GC /float(   self.raw_data_count * subdevide    )
			self.Q20_perc = 100*self.Q20/float(   self.raw_data_count * subdevide  )



			if self.trim:
				self.trim_data_count   = float(  self.trim_data_count )/subdevide
				self.trim_read_perc    = 100 * self.trim_read_count/float( self.raw_reads_count)
				self.trim_data_perc    = 100 * self.trim_data_count/self.raw_data_count

				output_rank  = [
				        self.raw_reads_count,
				        self.raw_data_count  ,
				        self.trim_read_count,
				        self.trim_data_count,
				        self.trim_read_perc ,
				        self.trim_data_perc,
				        self.filter_read_count,
				        self.filter_data_count,
				        self.filter_read_perc,
				        self.filter_data_perc,
				        self.N_perc,
				        self.GC_perc,
				        self.Q20_perc
				]
			else:
				output_rank  = [

				        self.raw_reads_count,
				        self.raw_data_count,
				        self.filter_read_count,
				        self.filter_data_count,
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

			raise  StopIteration
		
		output_hash[ 'pairname' ] = name[:-3]
		if self.trim:
			seq_re = re.search( '^N*([^N]\S{%s,}[^N])N*\n$'%( self.length-2 ) ,seq     )
			
			if  seq_re:
				[ trim_start,trim_stop ]     = seq_re.span(1)

				seq_trim      = seq[ trim_start:trim_stop  ]
				quality_trim  = quality[ trim_start : trim_stop  ]
				output_cont   = name+seq_trim+'\n'+define2+quality_trim+'\n'
				raw_cont            = output_cont
				seq                 = seq_trim
				quality             = quality_trim

				self.trim_data_count += len( seq_trim )
				self.trim_read_count += 1

				self.SUCCTRIM.write(  output_cont   )

				status = 'Trim'

				output_hash[ status ] = output_cont




			else:
				status = 'failTrim'

				output_hash[ status ] = raw_cont

				self.fail_trim.write(  raw_cont   )

				return (   output_hash     )


		if self.threshold_check(  quality   ):
			self.filter_read_count += 1
			self.filter_data_count += len( re.sub( '\s+','',  seq ) )
			#self.OUTPUT.write(  raw_cont  )
			status = 'Filter'
			output_hash[ status ] = raw_cont
		else:
			status = 'failFilter'
			output_hash[ status ] = raw_cont
			#self.FAIL_FILE.write(  raw_cont   )

		return output_hash


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
	
		## Split file into different group
	
		all_raw_file = glob.glob( os.path.abspath(path)+'/*.'+ options.appendix)
	
	
	
	
	
	
	
		trim = options.trim
	
		trim_length = options.trim_length
	
		quality = options.quality
		threshold = options.threshold
		return path , all_raw_file , quality , threshold , trim , trim_length
	except:
		print(  'The paramater you inpout is wrong !! please use -h to see the help!!!'  )
		sys.exit()
# run the proecss to qc
def SingleDbRun(  **paramater ):
	# the input variance is [  input_data , trim,  quality, trim_length ,threshold       ]
	QC_class = fastq_quality_class( **paramater  )
	for x in QC_class:
		pass
	#the output variance is ---> raw_reads( int,  1 ) , raw_data( float, Mb ), filter_reads(  int,1 ) , filter_data(  float,Mb  ), filter_perc(  float, 100*perc ), N_perc(  float, 100*perc ), GC_perc(  float, 100*perc ),Q20_perc(  float,100*perc ), FilterFilePath( char,  ), RawBoxGraph( fileName ),  RawDisGraph( fileName ) ,FilterBoxGraph( fileName ), FilterDisGraph( fileName )

	return QC_class.raw_reads_count , QC_class.raw_data_count, QC_class.filter_read_count ,QC_class.filter_data_count,QC_class.filter_data_perc, QC_class.N_perc, QC_class.GC_perc, QC_class.Q20_perc, QC_class.OUTPUT, QC_class.RawBoxGraph, QC_class.RawDisGraph, QC_class.FilterBoxGraph, QC_class.FilterDisGraph,QC_class.SUCCTRIM.name

def SingleEndRun( ):
	# to run qc to each file
	def quality_check(  each_f ):
		input_data = each_f
		data = {   }
		data .update( globals() )

		data.update( locals() )
		QC_class = fastq_quality_class( **data )
		for x in QC_class:
			pass


	other_paramater = locals()
	for each_f in all_raw_file:

		multiprocessing.Process( target =quality_check, args = ( each_f,  )   ).start()

if __name__ =='__main__':
	path,all_raw_file,quality,threshold,trim , trim_length = get_paramater()
	SingleEndRun(  )
