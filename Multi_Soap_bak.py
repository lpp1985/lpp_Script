#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/11
from lpp import *
import multiprocessing
import os,sys
from optparse import OptionParser 
import shutil, os
def check_path( path ):
	if os.path.exists(path):
		shutil.rmtree(path)
	os.mkdir( path )


def run(x):
	x=str(x)
	command = 'soapdenovo63mer all -s ./K'+x+'/config'+x+'.ini -K'+x+' -R -o ./K'+x+'/soap_k_'+x+' -p '+ str( cpu )+ ' >run'+x+'.log '
	os.system( command)
def assembly(   ):

	try:
		Kmer = eval( args[0] )
	except:
		print(  '''The paramater of Kmer must be a list,just like [ 1,2,3,4  ]  ''' )
		sys.exit()
	new = options.new
	if not Kmer:
		print( 'Without Kmer paramater, Programe could not execute!!'  )
		sys.exit()
	if new:
		for i in Kmer:
			check_path( 'K%s'%(i) )
	need_r = Kmer
	for i in need_r:
		CONFIG = open('./K%s/config%s.ini'%( i,i ),'w')
		CONFIG.write('''max_rd_len=100\n[LIB]\navg_ins=%s\nreverse_seq=0\nasm_flags=3\nrank=1\npair_num_cutoff=3\nrd_len_cutoff=%s\nmap_len=%s\nq1=%s\nq2=%s\n'''%( inszie,read_length,map_length,read1,read2 ))

		
	pool = multiprocessing.Pool(thread)
	pool.map( run, need_r)
	END = open( 'Assembly_status.output','w' )
	for e_f in glob.glob( '*.log' ):
		kmer = re.search( 'run(\d+)\.log',e_f ).group(1)
		
		RAW = open( e_f,'rU' )
		data = re.search( '\n(\d+\s+scaffolds from .+\n.+\n.+\n.+\n)' ,RAW.read() )
		if data:
			END.write( '%s\n>>>>>\n%s-----------------------------\n'%( kmer,data.group(1) )  )
if __name__=='__main__':
	usage = '''usage: python2.7 %prog [options] Kmer




	Kmer is a list of K value you want,e.g  [ 1, 2, 3, 4 ]'''
	parser = OptionParser(usage =usage ) 
	
	
	
	parser.add_option("-c", "--CPU", action="store", 
		              dest="cpu",
		              type='int', 
		              default = 15, 
		              help="CPU number for each thread")
	parser.add_option("-t", "--Thread", action="store", 
		              dest="thread", 
		              type='int', 
		              default = 3, 
		              help="CPU number for each thread")
	parser.add_option("-1", "--Reads1", action="store", 
		              dest="read1", 
		              help="Read1 location") 
	parser.add_option("-2", "--Reads2", action="store", 
		              dest="read2", 
		              help="Read2 location")
	
	parser.add_option("-n", "--new", action="store_true", 
		              dest="new", 
		              help="Do you want to build new directory to resore output,Appeared means true")
	parser.add_option("-i", "--inst", action="store", 
		              dest="inst", 
		              default = '400',
		              help="The insert size you want!!")
	parser.add_option("-m", "--map", action="store", 
		              dest="map", 
		              default = '30',
		              help="The Mapping length you want!!")
	parser.add_option("-r", "--readL", action="store", 
		              dest="readL", 
		              default = '60',
		              help="The Reads length you want!!")
	(options, args) = parser.parse_args() 
	cpu = options.cpu
	read1 = options.read1
	read2 = options.read2
	thread = options.thread
	inszie = options.inst
	map_length = options.map
	read_length = options.readL
	assembly()