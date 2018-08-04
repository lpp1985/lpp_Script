#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/7/13
from lpp import *
from optparse import OptionParser 
import sys
usage = '''usage: python2.7 %prog [options] 




for automatic statistic assembly result!!!'''
parser = OptionParser(usage =usage ) 
parser.add_option("-p", "--PATH", action="store", 
                  dest="path",
                  type='string', 
                  default = './', 
                  help="the assembly path")
parser.add_option("-o", "--OUTPUT", action="store", 
                  dest="output", 
                  type='string',  
                  help="The output file")
parser.add_option("-f", "--FILE", action="store", 
                  dest="file", 
                  type='string',  
                  help="The file name you want")


(options, args) = parser.parse_args() 

path = options.path

file_cate = re.escape(options.file)

OUTPUT = open( options.output ,'w'    )

log_hash = Ddict()
all_title = {}
for root,docr_list,file_list in os.walk(  path  ):
	query_name = root.split('/')[-1]
	for each_f in file_list:
		if re.search( '('+file_cate+')',each_f ):
			query_name = re.search( '(\d+)',each_f   ).group(1)
			LOG = open( root+'/'+each_f ,'rU'  )
			data = LOG.read()
			print( root )
			print( each_f )
			[(contig_totoal, contig_average , contig_n50, contig_n90       )] = re.findall( '''\d+ ctgs longer than 100, sum up (\d+)bp, with average length (\d+)
the longest is \d+bp, contig N50 is (\d+) bp,contig N90 is (\d+) bp''' ,data  )
			
			[ ( scaf_number , scaf_total,scaf_average,scaf_sing_number,scaf_sing_length , scaf_sing_average_length ,scaf_n50,scaf_n90 ) ] = re.findall( '''(\d+) scaffolds from \d+ contigs sum up (\d+)bp, with average length (\d+), \d+ gaps filled
(\d+) scaffolds&singleton sum up (\d+)bp, with average length (\d+)
the longest is \d+bp,scaffold N50 is (\d+) bp, scaffold N90 is (\d+) bp
Found 0 weak points in scaffolds''', data   )
			scaf_gap = re.search(  '(\d+) gaps overall' ,data  ).group(1)
			
			reads_in_gap_percentage = re.search( '''\((\S+)\)% reads in gaps''',data   ).group(1)
			reads_in_contig_percentage = re.search( '''\((\S+)\)% reads mapped to contigs''',data   ).group(1)
			reads_total = '%.2f'%( float(  reads_in_gap_percentage   ) + float(  reads_in_contig_percentage ) )+'%'
			reads_in_gap_percentage = reads_in_gap_percentage+'%'
			reads_in_contig_percentage = reads_in_contig_percentage+'%'
			variant = vars().copy()
			for each_key in variant:
				if  each_key.startswith('reads_')or each_key.startswith( 'scaf_'  ) or each_key.startswith( 'contig_'  ):
					exec(  '''log_hash[ query_name  ][ '%s' ] = \'%s\''''%( each_key, eval(each_key)  )  )
					print( '''log_hash[ query_name  ][ '%s' ] = \'%s\''''%( each_key, eval(each_key  ) )  )
					exec(   '''all_title[ '%s' ] = \'\''''%(  each_key  )  )
					print( query_name  , scaf_number)
		if each_f.endswith('.scafSeq'):
			query_name = re.search( '(\d+)',each_f   ).group(1)
			SCF = open( root+'/'+each_f   )
			seq_data = SCF.read()
			scaf_No_N_singl_number = str( len( re.findall( '([ATCG])' , seq_data   ) ) )
			scaf_N_singl_number = str(len( re.findall( '([N])' , seq_data   ) ))
			scaf_success = '%.4f'%( float( scaf_No_N_singl_number )*100/ int(log_hash[ query_name  ]['scaf_total'] ) )+'%'
			variant = vars().copy()
			for each_key in variant:
				if  each_key.startswith('scaf_N')or each_key.startswith( 'scaf_success'  ) :
					exec(  '''log_hash[ query_name  ][ '%s' ] = %s'''%( each_key, each_key  )  )
					exec(   '''all_title[ '%s' ] = \'\''''%(  each_key  )  )
OUTPUT.write(   'name'+'\t'+'\t'.join(  sorted( all_title  ) ) +'\n' )
for each_sample in sorted(  log_hash  ):
	OUTPUT.write( each_sample  )
	for each_data in sorted( log_hash[ each_sample ] ):
		OUTPUT.write( '\t'+log_hash[ each_sample ][each_data]  )
		
	OUTPUT.write( '\n' )
