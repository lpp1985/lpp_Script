#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: Filter bad quality reads for denovo sequencing
# Created: 2011/5/9
'''usage quality_stats.py [path]  [appendix]'''
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
parser.add_option("-f", "--FILE", action="store", 
                  dest="file", 
                  help="appendix name represent sequencing data")
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
for prefix in  name_hash:
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

