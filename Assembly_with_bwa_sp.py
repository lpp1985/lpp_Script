#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/9
from lpp import *
from optparse import OptionParser
from Multi_Soap import *
import Multi_Soap
'# you could type [  SCRIPT_NAME  ] -h to see the documentation !!!!'
usage='''usage: python %prog [options] 

It can assembly reads from the infered insert size of BWA'''

parser = OptionParser(usage =usage ) 

parser.add_option("-f", "--FILE", action="store", 
                  dest="bwa",
                  default = 'log',
                  type='string',  
                  help="Input BWA SAMPE LOG")
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

parser.add_option("-m", "--map", action="store", 
                  dest="map", 
                  default = '30',
                  help="The Mapping length you want!!")
parser.add_option("-r", "--readL", action="store", 
                  dest="readL", 
                  default = '60',
                  help="The Reads length you want!!")
(options, args) = parser.parse_args() 

'''path is the PATH of bwa sampe log which print the inset size '''
'''appendix is the appendix name of SAMPE LOG'''

LOG_FILE = open( os.path.abspath( options.bwa ) )
def run(x):
	x=str(x)
	command = 'soapdenovo63mer all -s ./K'+x+'/config'+x+'.ini -K'+x+' -R -o ./K'+x+'/soap_k_'+x+' -p '+ str( cpu )+ ' >run'+x+'.log '
	os.system( command)
inszie = re.search( 'inferred maximum insert size: (\d+) ',  LOG_FILE.read()  ).group(1)
''' insert_size is the size of data '''
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

if __name__ =='__main__':
	cpu = options.cpu
	read1 = options.read1
	read2 = options.read2
	thread = options.thread
	map_length = options.map
	read_length = options.readL
	assembly()