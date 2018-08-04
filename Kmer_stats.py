#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/9
'''usage quality_stats.py [path]  [appendix]'''
from lpp import *
'''usage: quality_stats.py [reads] [ prefix ]'''
import multiprocessing
#def fastx_check(  reads,end_path ):

	#end_title = end_path + reads

	#os.system(  '''fastq_quality_filter -q 20 -p 75 -i %s -o %s.filter  >%s.log
#'''%( reads,end_title,end_title   )   )
	#os.system(  '''fastx_quality_stats -i %s.filter > %s.qualstats'''%(  end_title,end_title+'.filter'  )   )
	#os.system(  '''fastx_quality_boxplot.R %s.qualstats %s.qualstats.png'''%(  end_title+'.filter',end_title+'.filter'  )   )
	#os.system('''fastx_nucleotide_distribution.R  %s.qualstats %s.nucdistr.png'''%(  end_title+'.filter',end_title+'.filter'   ))
	#os.system(  '''fastx_quality_stats -i %s > %s.qualstats'''%(  reads,end_title  )   )
	#os.system(  '''fastx_quality_boxplot.R %s.qualstats %s.qualstats.png'''%(  end_title,end_title  )   )
	#os.system('''fastx_nucleotide_distribution.R  %s.qualstats %s.nucdistr.png'''%(  end_title,end_title   ))

def corrected_check( reads,end_path ):
	end_title = end_path + reads
	os.system(  '''fastx_quality_stats -i %s > %s.qualstats'''%(  end_title,end_title  )   )
	os.system(  '''fastx_quality_boxplot.R %s.qualstats %s.qualstats.png'''%(  end_title,end_title  )   )
	os.system('''fastx_nucleotide_distribution.R  %s.qualstats %s.nucdistr.png'''%(  end_title,end_title   ))



input_path = sys.argv[1]
file_appendix = sys.argv[2]
all_f = glob.glob( input_path+'/*.'+ file_appendix )
name_hash = Ddict()
for each_f in all_f:
	f_name = os.path.split(each_f)[-1].split('.')[0]
	name_hash[ f_name ] [ each_f ] = ''
for prefix in  name_hash:
	end_path = input_path+'/Kmer_freq/'+prefix+'/'
	if not os.path.isdir(  end_path ):
		os.makedirs( end_path )
	NEED_LIST = open(  end_path+prefix+'.list','w' )

	for each_reads in name_hash[ prefix ] :
		NEED_LIST.write( input_path+'/'+ each_reads +'\n'  )
		fastx_check( each_reads, end_path )



	os.system( 'KmerFreq -i %s -o ./%s  >kmer_fre_error'%( NEED_LIST.name,end_path+prefix   )   )
	os.system( '''Corrector -i %s -r %s%s.freq -t 20'''%( NEED_LIST.name ,end_path,prefix   )  )
	for each_reads in name_hash[ prefix ] :
		if 'read1' in each_reads:
			read1 = each_reads
		else:
			read2 = each_reads

	os.system(  '''extract_pair.py  %s%s %s%s %s %s %s'''%( input_path,read1+'.corr',input_path,read2+'.corr' ,end_path + prefix + '.corrected' ,read1+'.filter' ,read2+'.filter' )                                            )
	all_cooreted = glob.glob( K_merdir+'*.fastq'   )
	for each_f in all_cooreted:

		corrected_check( each_f,'./' )

