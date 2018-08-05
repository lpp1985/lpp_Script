#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
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



def fastx_filter_check(  reads,end_path ):
	
	end_title = reads
	
	os.system(  '''fastx_quality_stats -i %s > %s.qualstats'''%(  reads,end_title  )   )
	
	os.system(  '''fastx_quality_boxplot.R %s.qualstats %s.qualstats.png'''%(  end_title,end_title  )   )

	
	
	os.system('''fastx_nucleotide_distribution.R  %s.qualstats %s.nucdistr.png'''%(  end_title,end_title   ))
	
	if 'read1' in reads:
		RAW = fastq_check(  open( reads,'rU' ) )
		length_reads = len( RAW.next()[1][:-1]  )
		STATS = open( '%s.stats'%( end_title ),'w' )
		STATS.write(  os.popen(  '''aa=`wc -l %s| grep -Po '^\d+'` &&echo "scale=3;$aa/2*%s/1024/1024"|bc  && echo "$aa/2"|bc  '''%( end_title , length_reads  )  ).read()  )
	
	
	
	
def fastx_check(  reads,end_path ):

	end_title = end_path + reads
	
	os.system(  '''fastx_quality_stats -i %s > %s.qualstats'''%(  reads,end_title  )   )
	
	os.system(  '''fastx_quality_boxplot.R %s.qualstats %s.qualstats.png'''%(  end_title,end_title  )   )

	
	
	os.system('''fastx_nucleotide_distribution.R  %s.qualstats %s.nucdistr.png'''%(  end_title,end_title   ))
	
	os.system(  '''fastq_quality_filter -q 20 -p 80 -i %s -o %s.other  >%s.log'''%( reads,end_title,end_title   )   )
	
	
	if 'read1' in reads:
		RAW = fastq_check(  open( reads,'rU' ) )
		length_reads = len( RAW.next()[1][:-1]  )
		STATS = open( '%s.stats'%( end_title ),'w' )
		STATS.write(  os.popen(  '''aa=`wc -l %s| grep -Po '^\d+'` &&echo "scale=3;$aa/2*%s/1024/1024"|bc  && echo "$aa/2"|bc  '''%( reads , length_reads  )  ).read()  )
	
	#os.system(  '''fastx_quality_stats -i %s.filter > %s.qualstats'''%(  end_title,end_title+'.filter'  )   )
	
	
	#os.system(  '''fastx_quality_boxplot.R %s.qualstats %s.qualstats.png'''%(  end_title+'.filter',end_title+'.filter'  )   )
	
	
	#os.system('''fastx_nucleotide_distribution.R  %s.qualstats %s.nucdistr.png'''%(  end_title+'.filter',end_title+'.filter'   ))
	
	

input_path = options.path
file_appendix = options.file
all_f = glob.glob( input_path+'/*.'+ file_appendix )
name_hash = Ddict()
for each_f in all_f:
	f_name = os.path.split(each_f)[-1].split('.')[0]
	name_hash[ f_name ] [ each_f ] = ''

for prefix in  name_hash:
	end_path = input_path+'/Filter/'+prefix+'/'
	if not os.path.isdir(  end_path ):
		os.makedirs( end_path )

	for each_reads in name_hash[ prefix ] :
		
		fastx_check( each_reads, end_path )
	if len( name_hash[ prefix ]     )>1:
		read1 =[  x for x in     name_hash[ prefix ]  if 'read1' in x    ] [ 0 ] 
		read2 =[  x for x in     name_hash[ prefix ]  if 'read2' in x    ] [ 0 ]
		print(  read1 )
	os.system( ''' check_pair.py -1 %s -2 %s -o %s    '''%( end_path+read1+'.other',   end_path+read2+'.other', end_path+re.search('(\S+)\.read\d+\.fastq',read1).group(1)+'.filter'     )      )
	for key1 in glob.glob(  end_path+'*.filter.read?.fastq'    ):
		fastx_filter_check( key1,end_path   )
	
