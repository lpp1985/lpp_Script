#!/usr/bin/env python
#coding:utf-8
# Author:  LPP
# Purpose:
# Created: 2011/11/7
from lpp import *
from datetime import datetime
from time import strftime
os.system(  'source ~lipp/.bashrc'  )
import multiprocessing
from  optparse import OptionParser
import subprocess
import os,re
from quality_stats_Pairend_Distri import get_paramater,PairEndRun,get_sample,fastq_quality_pair
SYS = os.getenv( 'SYS'  )
import sys
#sys.path.append( SYS )
#from Project.models import *

def bowtie_check(  each_group  ):
	'''do bowtie and output the list to build expression matrix'''
	sample_name = os.path.split(  all_file_group[ each_group  ][0]  )[-1].split('.')[0]
	if len(  all_file_group[ each_group  ] ) >1:

		[read1file ,read2file] =  sorted(
			    all_file_group[ each_group  ],
			    key = lambda x:
			    int(
			            re.search( 'read(\d+)' ,x ).group(1)
			          )
			    )
		[  read1name,read2name   ]  = map(  lambda x : os.path.split(x)[-1],   [read1file ,read2file]       )

		bowtie_subprocess = subprocess.Popen(
			    [
			            'bowtie',
			            '-X',       '1000',
			            '-I',            '50',
			            '-c' ,'%s'%( enddata ),
			            '--chunkmbs', '1000',
			            '-a',
			            '-v', '3',
			            '-M' ,'1',
			            '--best',
			            '--strata',
			            '-q',
			            '-1' ,read1file,
			            '-2',read2file ,
			            '--trim3', '16',
			            '-e', '120',
			            '--phred64-quals',
			            '-p', '%s'% (  thread  ),
			            '--al', sample_name+'.Mappable',
			            '--un', sample_name+'.Un',
			            sample_name+'.list',


			    ],
			    stdout=subprocess.PIPE,
			    stderr = subprocess.PIPE,
		)
	else:
		read_file = all_file_group[ each_group  ][0]
		bowtie_subprocess = subprocess.Popen(
				        [
				                'bowtie',
				                '-c' ,'%s'%( enddata ),
				                '--chunkmbs', '1000',
				                '-a',
				                '-v', '3',
				                '-M' ,'1',
				                '--best',
				                '--strata',
				                '-q',
		                        read_file,
				                '--trim3', '16',
				                '-e', '120',
				                '--phred64-quals',
				                '-p', '%s'% (  thread  ),
				                '--al', sample_name+'.Mappable',
				                '--un', sample_name+'.Un',
				                sample_name+'.list',


				        ],
				        stdout=subprocess.PIPE,
				        stderr = subprocess.PIPE,
				)

	stdout,stderr = bowtie_subprocess.communicate()
	subprocess.call(
	        [
	                'analysis_bowtie.py',
	                sample_name+'.list',
	                sample_name+'.count'
	         ],
	        )
	print( '''-------------%s 's mapping status-----------------------'''%(  each_group  )  )
	print(  stderr+'\n'  )
	print(   '''--------------------------------------------------------------------'''  )



def gicl_assembly( inputdata   ):
	inputdata = os.path.abspath(  inputdata )







	os.system( ' '.join(
		        [
		                'cd' ,os.path.dirname(  inputdata  ),
		                '&&',
		                'gicl',
		                '%s'%(  inputdata ),
		                '-c', '10',
		                '-n', '1000',
		                '-p' ,'95',
		                '-v',' 5',
		                '-l', '50',
		                ],
		)
		                )

	ace_file = inputdata+'_self.ace'
	singletonSeq = inputdata+'.Singleton'
	assemblyend = inputdata+'.end'
	subprocess.call(   [ 'get_seq.py',
	                     '-s', '%s'%(  inputdata ),
	                     '-l',inputdata+'_self.singletons' ,
	                     '-n', '1',
	                     '-o' ,singletonSeq,
	                     ],   )
	subprocess.call(   [ 'gicl_combine.py',
	                     '-a', '%s'%(  ace_file ),
	                     '-t','1',
	                     '-o', assemblyend,
	                     '-s' ,singletonSeq,
	                     ],   )
	return os.path.getsize( ace_file), assemblyend




#check and build new folder

def SignificantExpress( inputpath, outputpath   ):
	check_path( outputpath  )

	stderr,stdout= subprocess.Popen(
	        [
	                'Expression_different_R_static.py',
	                '-o',outputpath,
	                '-a','matrix',
	                '-s',path+'stats/',
	                '-i',inputpath,
	                '-t','0.5',
	                ],
	        stdout=subprocess.PIPE,
	        stderr = subprocess.PIPE

	        ).communicate()
	print( stderr )
def Path_checking( path ):
	path = os.path.abspath( path )
	if os.path.exists(  path ):
		shutil.rmtree( path )
	os.makedirs(  path  )
	return path+'/'
def getpara(  ):
	'''get the paramater to run the programe'''
	try:
		'# you could type [  SCRIPT_NAME  ] -h to see the help log !!!!'
		usage='''usage: python %prog [options]

	    It can automaticaly analysis transcriptome sequence'''
		parser = OptionParser(usage =usage )
		parser.add_option("-p", "--Path", action="store",
			          dest="path",
			          type='string',
			          help="Input path")
		#parser.add_option("-P", "--Project", action="store",
				                  #dest="project",
				                  #type='string',
				                  #help="Project")
		parser.add_option("-a", "--Appendix", action="store",
			          dest="appendix",
			          type='string',
			          default='fastq',
			          help="Input FASTQ")

		parser.add_option("-f", "--FINSH_Appendix", action="store",
				                  dest="finish",
				                  type='string',
				                  default='pair*',
				                  help="Input FASTQ")





		#parser.add_option("-t", "--Trim", action="store_true", default= False,
			          #dest="trim",
			          #help="Do you want to trim the fastq file( trim the N in head or in tail  )")


		#parser.add_option("-l", "--Trim_length", action="store",
			          #dest="trim_length",
			          #type= "int",
			          #default=60,
			          #help="The length of threshold  to  filter after trim, THIS PARAMATER OUGHT TO BE INGORANCE IN [-t] IS MISSING!!!")


		#parser.add_option("-q", "--Quality", action="store", type = "int",
			          #dest="quality",
			          #default=20,
			          #help="Minimum quality score to keep.")


		#parser.add_option("-r", "--Threshold", action="store", type = "float",
			          #dest="threshold",
			          #default = 0.2,
			          #help='''The Threashold you want to filter, if this paramater's value locate in the area of (0,1), It represent the maxium percentage the low quality accounts for ,
			      #if Integer it represents the maxium number of low quality base in the whole reads to trim''')
		parser.add_option("-P", "--Process", action="store",
			          dest="process",
			          type='int',
			          help="process Number")
		parser.add_option("-T", "--Thread", action="store",
			          dest="thread",
			          type='int',
			          default=20,
			          help="Thread number")
		parser.add_option("-C", "--Cat", action="store",
			          dest="catlength",
			          type='int',
			          default=200,
			          help="The threshold number to filter cat_all")
		#parser.add_option("-b", "--Block", action="store",
				                                  #dest="blocks",
				                                  #type='int',
				                                  #default=50000,
				                                  #help="each block size to map reduce")



		#parser.add_option("-x", "--Cache", action="store",
                                                                  #dest="cache",
                                                                  #type='int',
                                                                  #default=10,
                                                                  #help="The cache size of to Output")

		(options, args) = parser.parse_args()

		path = options.path

		## Split file into different group

		all_raw_file = glob.glob( os.path.abspath(path)+'/*.'+ options.appendix)
		all_finished_file=[]
		for a,b,c in os.walk(  os.path.abspath(path) ):
			for each_f in c:
				if re.findall( '(\.'+options.finish+'*)',each_f  ):
					all_finished_file.append( a+'/'+each_f    )
			
		#all_finished_file = glob.glob( os.path.abspath(path)+'/*.'+ options.finish )
		







		#trim = options.trim

		#trim_length = options.trim_length

		#quality = options.quality
		#threshold = options.threshold
		process = options.process
		thread = options.thread
		catlength = options.catlength
		#cache_size = options.cache
		#block_size = options.blocks
		return path ,all_finished_file, all_raw_file, process, thread,catlength
	except:

		print(  'The paramater you inpout is wrong !! please use -h to see the help!!!'  )
		sys.exit()



def QC( each_group , QUEUE  ):
	'''do QC and return the filtered file name'''
	if len(  all_file_group[ each_group ] )==2:

		print(   '%s start to QC'%( each_group )  )
		data = globals()
		files = all_file_group[ each_group ]
		data.update( locals()  )
		qc_instance = fastq_quality_pair( **data )
		qc_instance.Qc()
		QUEUE.put( [qc_instance.FILTERPAIR1.name, qc_instance.FILTERPAIR2.name] )
	#return qc_instance.FILTERPAIR1.name, qc_instance.FILTERPAIR2.name

def Total_Run(  each_group ):
	QUEUE = multiprocessing.Queue(  )
	#Qc_process = multiprocessing.Process(  target= QC,args=[   each_group   , QUEUE]  )
	#Qc_process.daemon = False
	#Qc_process.start()
	#Qc_process.join()
	#if QUEUE.qsize():
		#Filter_read1, Filter_read2 = QUEUE.get()
	#else:
		#raise IOError ,'  %s Qc process has error!! Please Check!! '
	#print(  '%s QC ok! Start to Soap test' %( each_group  ) )
	if len(  all_finished_group[ each_group ] ) >1:
		[Filter_read1, Filter_read2 ]= sorted(
		        all_finished_group[ each_group ] ,
		        key= lambda x: int( re.search( 'pair(\d+)',x  ).group(1)  )
		                                       )
		Step2_prePath = Path_checking( '%s../Step2_pre/%s/'%(path, each_group )  )



		TestSoap_process = multiprocessing.Process(
			target= Soap,
			                                     kwargs=
			                                     {
			                                             'read1file' :Filter_read1,
			                                             'read2file' : Filter_read2,
			                                             'outputpath':Step2_prePath,
			                                             'Queue' : QUEUE
			                                     }
		)
		TestSoap_process.daemon = False
		TestSoap_process.start()
		TestSoap_process.join()
		data = QUEUE.get()
		if 'IOError' in data:
			raise IOError , data[-1]
		else:
			contigRef = data
		print( '%s Soap test ok! try to BWA'%( each_group ) )
		Step3_BWAinferPath = Path_checking( '%s../Step3_BWAInfer/%s/'%( path,each_group )  )
		BWA_process = multiprocessing.Process(  target=BWA_Run,
			                                kwargs=  {
			                                        'ref'                 :     contigRef,

			                                        'inputpath'    :      os.path.dirname( Filter_read1    ),
			                                        'outputpath' :       Step3_BWAinferPath,
			                                        'Queue'         :        QUEUE
			                                        }
			                                )
		BWA_process.start()
		BWA_process.join()
		data = QUEUE.get()
		if 'IOError' in data:
			raise IOError , data[-1]
		else:
			insertsize = data
		print(  '%s bwa ok! insertsize = %s , Try to Multi Soap!!!'%( each_group,insertsize ) )
		Step4_MultiSoap_Path = Path_checking( '%s../Step4_SoapMulti/%s/'%( path , each_group )  )
		MultiSoap_process = multiprocessing.Process(
			target= Soap,
			                                     kwargs=
			                                     {
			                                             'read1file' :Filter_read1,
			                                             'read2file' : Filter_read2,
			                                             'outputpath':Step4_MultiSoap_Path,
			                                             'Queue' : QUEUE    ,
			                                             'insertsize':insertsize,
			                                             'kmer':[ 21,25,29,33,37,41,45,49,53,57  ]
			                                     }
		)
		MultiSoap_process.start()
		MultiSoap_process.join()
		print( '%s Soap ok!! Finish!!' %( each_group ))

def Soap(  read1file,read2file ,Queue, outputpath ,insertsize=400,kmer=['21',''] ):
	soap_process = subprocess.Popen([  'Multi_Soap.py',
	                                   '-c','%s'%( thread )     ,
	                                   '-t','%s'%( process )    ,
	                                   '-1','%s'%( read1file )  ,
	                                   '-2','%s'%( read2file )  ,
	                                   '-n',
	                                   '-i','%s'%( insertsize  ),
	                                   '-o','%s'%( outputpath  ),
	                                   ','.join(map( str ,kmer  )  ),
	                                   ]
	                                ,
	                                stdout= subprocess.PIPE,
	                                stderr= subprocess.PIPE
	                                )
	stdout,stderr = soap_process.communicate()
	#if stderr:
		#Queue.put( ['IOError',stderr] )
		#raise IOError,stderr
	for a,b,c in os.walk( outputpath ):
		for each_f in c:
			if each_f.endswith( '.contig' ):
				Queue.put( a+'/'+each_f )

def BWA_Run( ref,inputpath,outputpath , Queue   ):
	print( ' '.join(   [
	                'BWA_Mapping.py',
	                '-r',         ref        ,
	                '-t','%s'%( thread )     ,
	                '-p','%s'%( process )    ,
	                '-o',    outputpath      ,
	                '-i',  inputpath         ,
	                ]  )  )
	bwa_process = subprocess.Popen(
	        [
	                'BWA_Mapping.py',
	                '-r',         ref        ,
	                '-t','%s'%( thread )     ,
	                '-p','%s'%( process )    ,
	                '-o',    outputpath      ,
	                '-i',  inputpath         ,
	                ],
	        stdout= subprocess.PIPE,
	        stderr= subprocess.PIPE
	        )
	stdout,stderr = bwa_process.communicate()
	if stderr:
		Queue.put( stderr )
		sys.exit()
	else:
		for a,b,c in os.walk( outputpath ):
			for each_f in c :
				if each_f.endswith('_PE.log'):
					all_data = open(  a+'/'+each_f ,'rU'  ).read()
					infer_size = re.search(  'percentile: \(\d+, (\d+), \d+\)',all_data ).group(1)
					Queue.put(  infer_size )
					return infer_size


if __name__=='__main__':
	# get the paramater
	path ,all_finished_file, all_raw_file, process, thread,catlength= getpara()

	path = os.path.abspath( path )
	if path[-1]!='/':
		path+='/'

	all_file_group = {}
	all_finished_group = {}
	#seperate the sample to different group
	for each_f in all_raw_file:
		name = get_sample( each_f )
		name_suffix = re.escape( name )
		all_file_group[ name] = [  x   for x in all_raw_file if re.search( '^('+name_suffix+')',  os.path.split(os.path.abspath(x)  )[-1]  )     ]

	for each_f in all_finished_file:
		name = get_sample( each_f )
		name_suffix = re.escape( name )
		all_finished_group[ name] = [  x   for x in all_finished_file if re.search( '^('+name_suffix+')',  os.path.split(os.path.abspath(x)  )[-1]  )     ]
	print(  all_file_group )
	print(  all_finished_group )
	# get the sample

	all_raw_assembly_process = []
	for  each_group in all_finished_group:

		Total_process = multiprocessing.Process(  target=Total_Run,args=[  each_group  ,]  )
		Total_process.daemon = False
		Total_process.start()
		all_raw_assembly_process.append( Total_process )
	for each_raw_process in all_raw_assembly_process:
		each_raw_process.join()
	print('''######################################################
######################################################
######################################################
	''')
	print( 'All the Multiple Ker Assembly is Finished ,Start to cluster and Assembly'  )

	Step5_path = Path_checking( '%s../Step5_GICL/'%(path)  )

	os.system( 'cat_all_fasta.py -o %sall_contig -a contig -i %s../Step4_SoapMulti/'%(  Step5_path  ,   path   )  )
	needgiclfile ='%sall_contig.fasta'%( Step5_path )






	times = 1
	ace_size = 1
	while ace_size:
		containerPath = Step5_path+'%s/'%( times )
		if  os.path.exists(  containerPath  ):
			shutil.rmtree( containerPath )
		os.makedirs( containerPath )

		inputdata = containerPath+str(  times  )
		os.system(  'ln -s %s  %s'%(  needgiclfile  ,   inputdata      )       )
		ace_size, needgiclfile  = gicl_assembly(  inputdata   )

		times+=1
		enddata = Step5_path +'FinalEnd.fasta'
	END = open( enddata,'w'  )
	for t,s in fasta_check( open(  needgiclfile ,'rU')   ):
		s1 = re.sub( '\s+','',s )
		if len( s1) >= catlength:
			END.write( t+s )


	print(  'Assembly finished!!!' )
	# now to bowtie
	print( 'Start to bowtie and analysis expression profile' )
	Step6_bowtie = Path_checking( '%s../Step6_Bowtie_matrix/'%(path)  )
	# build bowtie index

	stdout, stderr  = subprocess.Popen(

	        [
	                'bowtie-build','-f',
	                enddata,
	                enddata,
	                '2>&1',
	                '1>/dev/null',
	                '3>&1',
	                ],
	        stdout=subprocess.PIPE,
	        stderr= subprocess.PIPE

	).communicate()

	all_bowtie_process = []
	os.chdir( Step6_bowtie )
	for each_group in all_file_group:

		bowtie_process = multiprocessing.Process( target=bowtie_check, args=[ each_group,  ]   )
		bowtie_process.daemon = False
		bowtie_process.start()
		all_bowtie_process.append(  bowtie_process )
	for each_process in all_bowtie_process:
		each_process.join()
	# build expression matrix
	os.system(
	                'expression_matrix_build.py   *.count'

	                 )
	# do significant test
	Step7_Deseq = '%s../Step7_Deseq/'%(path)
	SignificantExpress( Step6_bowtie , Step7_Deseq  )