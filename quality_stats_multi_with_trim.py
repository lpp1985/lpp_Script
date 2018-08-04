#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/8/29
from lpp import *
import multiprocessing
from optparse import OptionParser
usage='''usage: python %prog [options] 

It can stats quality and filter bad quality reads'''
parser = OptionParser(usage =usage ) 
parser.add_option("-p", "--Path", action="store", 
                  dest="path",
                  type='string',  
                  help="Input path")
parser.add_option("-a", "--appendix", action="store", 
                  dest="file", 
                  help="appendix name represent sequencing data")
#parser.add_option("-i", "--Input", action="store", 
                  #dest="input",
                  #type='string',  
                  #help="Input FASTQ")

#parser.add_option("-o", "--OUTPUT", action="store", 
                  #dest="output", 
                  #help="The fastq file successfully filter the quality")


#parser.add_option("-f", "--fail", action="store", 
                  #dest="fail_file", 
                  #help="The fastq file fail to  filter the quality")


parser.add_option("-t", "--Trim", action="store_true", default= False,
                  dest="trim", 
                  help="Do you want to trim the fastq file( trim the N in head or in tail  )")

#parser.add_option("-j", "--Fail_Trim", action="store", 
                  #dest="fail_trim", 
                  #help="The reads fail to  satified the trim_length_threshold, THIS PARAMATER OUGHT TO BE INGORANCE IN [-t] IS MISSING!!!")

#parser.add_option("-s", "--Success_Trim", action="store", 
                  #dest="succ_trim", 
                  #help="The reads successfully   satified the trim_length_threshold, THIS PARAMATER OUGHT TO BE INGORANCE IN [-t] IS MISSING!!!")
parser.add_option("-l", "--Trim_length", action="store", 
                  dest="trim_length", 
                  type= "int",
                  help="The length of threshold  to  filter after trim, THIS PARAMATER OUGHT TO BE INGORANCE IN [-t] IS MISSING!!!")


parser.add_option("-q", "--Quality", action="store", type = "int",
                  dest="quality", 
                  help="Minimum quality score to keep.")


parser.add_option("-r", "--Threshold", action="store", type = "float",
                  dest="threshold", 
                  help='''The Threashold you want to filter, if this paramater's value locate in the area of (0,1), It represent the maxium percentage the low quality accounts for , 
                  if Integer it represents the maxium number of low quality base in the whole reads to trim''')
print(  sys.argv[1:]  )
(options, args) = parser.parse_args() 


def fastx_check(  reads,end_path ):

	end_title = end_path + reads

	os.system(  '''fastq_quality_filter -q 20 -p 80 -i %s -o %s.filter  >%s.log
'''%( reads,end_title,end_title   )   )
	os.system(  '''fastx_quality_stats -i %s.filter > %s.qualstats'''%(  end_title,end_title+'.filter'  )   )
	os.system(  '''fastx_quality_boxplot.R %s.qualstats %s.qualstats.png'''%(  end_title+'.filter',end_title+'.filter'  )   )
	RAW = fastq_check(  open( reads,'rU' ) )
	length_reads = len( RAW.next()[1][:-1]  )
	if 'read1' in reads:
		STATS = open( '%s.stats'%( end_title ),'w' )
		STATS.write(  os.popen(  '''aa=`wc -l %s| grep -Po '^\d+'` &&echo "scale=3;$aa/2*%s/1024/1024"|bc  && echo "$aa/2"|bc  '''%( reads , length_reads  )  ).read()  )
		#print( '''aa=`wc -l %s| grep -Po '^\d+'` &&echo "scale=3;$aa/2*%s/1024/1024/1024"|bc  && echo "$aa/2"|bc  >%s.stats'''%( end_title+'.filter' ,length_reads, end_title+'.filter'  ))
		FIL_STATS = open( '%s.stats'%( end_title+'.filter' )  ,'w'  )
		FIL_STATS.write(  os.popen(  '''aa=`wc -l %s| grep -Po '^\d+'` &&echo "scale=2;$aa/2*%s/1024/1024"|bc  && echo "$aa/2"|bc  '''%( end_title+'.filter' ,length_reads  )  ).read()  ) 
	os.system('''fastx_nucleotide_distribution.R  %s.qualstats %s.nucdistr.png'''%(  end_title+'.filter', end_title+'.filter'  ))
	os.system(  '''fastx_quality_stats -i %s > %s.qualstats'''%(  reads,end_title  )   )
	os.system(  '''fastx_quality_boxplot.R %s.qualstats %s.qualstats.png'''%(  end_title,end_title  )   )
	os.system('''fastx_nucleotide_distribution.R  %s.qualstats %s.nucdistr.png'''%(  end_title,end_title   ))

def corrected_check( reads,end_path ):
	end_title = end_path + reads
	RAW = fastq_check(  open( reads,'rU' ) )
	length_reads = len( RAW.next()[1][:-1]  )
	os.system(  '''fastx_quality_stats -i %s > %s.qualstats'''%(  end_title,end_title  )   )
	os.system(  '''fastx_quality_boxplot.R %s.qualstats %s.qualstats.png'''%(  end_title,end_title  )   )
	os.system('''fastx_nucleotide_distribution.R  %s.qualstats %s.nucdistr.png'''%(  end_title,end_title   ))
	if 'pair1' in reads:
		K_STATS = open( end_title+'.stats'      ,'w'  )
		
		K_STATS.write(  os.popen(  '''aa=`wc -l %s| grep -Po '^\d+'` &&echo "scale=2;$aa/2*%s/1024/1024"|bc  && echo "$aa/2"|bc'''%( reads ,length_reads )  ) .read()  )
	
	
input_path = options.path
file_appendix = options.file
all_f = glob.glob( input_path+'/*.'+ file_appendix )
name_hash = Ddict()
for each_f in all_f:
	f_name = os.path.split(each_f)[-1].split('.')[0]
	name_hash[ f_name ] [ each_f ] = ''

if not os.path.isdir(  './Kmer_freq/' ):
	os.makedirs( './Kmer_freq/' )

#for prefix in  name_hash:
	#end_path = input_path+'/Filter/'+prefix+'/'
	#if not os.path.isdir(  end_path ):
		#os.makedirs( end_path )
	#NEED_LIST = open(  end_path+prefix+'.list','w' )

	#for each_reads in name_hash[ prefix ] :
		#NEED_LIST.write( end_path+ each_reads +'.filter\n'  )
		#fastx_check( each_reads, end_path )

	#K_merdir =  './Kmer_freq/'+prefix +'/'
	#os.mkdir(  K_merdir )	

	#os.system( 'KmerFreq -i %s -o ./%s  >kmer_fre_error'%( NEED_LIST.name,K_merdir+prefix   )   )
	#os.system( '''Corrector -i %s -r %s%s.freq -t 20'''%( NEED_LIST.name ,K_merdir,prefix   )  )
	#for each_reads in name_hash[ prefix ] :
		#if 'read1' in each_reads:
			#read1 = each_reads
		#else:
			#read2 = each_reads
	
	#os.system(  '''extract_pair.py  %s%s %s%s %s %s%s %s%s'''%( end_path,read1+'.filter.corr',end_path,read2+'.filter.corr' ,K_merdir+prefix+'.corrected' ,end_path,read1+'.filter' ,end_path,read2+'.filter' )                                            )
	#all_cooreted = glob.glob( K_merdir+'*.fastq'   )
	#for each_f in all_cooreted:

		#corrected_check( each_f,'./' )
def extract_quality ( prefix  ):
	end_path = input_path+'/Filter/'+prefix+'/'
	if not os.path.isdir(  end_path ):
		os.makedirs( end_path )
	NEED_LIST = open(  end_path+prefix+'.list','w' )

	for each_reads in name_hash[ prefix ] :
		NEED_LIST.write( end_path+ each_reads +'.filter\n'  )
		fastx_check( each_reads, end_path )

	K_merdir =  './Kmer_freq/'+prefix +'/'
	os.mkdir(  K_merdir )	

	os.system( 'KmerFreq -i %s -o ./%s  >kmer_fre_error'%( NEED_LIST.name,K_merdir+prefix   )   )
	os.system( '''Corrector -i %s -r %s%s.freq -t 20'''%( NEED_LIST.name ,K_merdir,prefix   )  )
	for each_reads in name_hash[ prefix ] :
		if 'read1' in each_reads:
			read1 = each_reads
		else:
			read2 = each_reads
	
	os.system(  '''extract_pair.py  %s%s %s%s %s %s%s %s%s'''%( end_path,read1+'.filter.corr',end_path,read2+'.filter.corr' ,K_merdir+prefix+'.corrected' ,end_path,read1+'.filter' ,end_path,read2+'.filter' )                                            )
	all_cooreted = glob.glob( K_merdir+'*.fastq'   )
	for each_f in all_cooreted:

		corrected_check( each_f,'./' )
pool = multiprocessing.Pool( processes=20 ) 
pool.map( extract_quality,  name_hash )