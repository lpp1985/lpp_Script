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
parser.add_option("-3", "--Reads3", action="store", 
                  dest="read3", 
                  help="Read3 location") 
parser.add_option("-4", "--Reads4", action="store", 
                  dest="read4", 
                  help="Read4 location")

parser.add_option("-n", "--new", action="store_true", 
                  dest="new", 
                  help="Do you want to build new directory to resore output,Appeared means true")
parser.add_option("-i", "--inst1", action="store", 
                  dest="inst1", 
                  default = '400',
                  help="The insert size you want!!")
parser.add_option("-I", "--inst2", action="store", 
                  dest="inst2", 
                  default = '400',
                  help="The dataset2's insert size you want!!")
parser.add_option("-s", "--Rversed_seq1", action="store", 
                  dest="s1", 
                  default = '0',
                  help="The dataset1's reversed_seq you want!![ Pair-end: 0, Mate-pair : 1   ]")
parser.add_option("-S", "--Rversed_seq2", action="store", 
                  dest="s2", 
                  default = '0',
                  help="The dataset2's reversed_seq you want!![ Pair-end: 0, Mate-pair : 1   ]")
parser.add_option("-m", "--map1", action="store", 
                  dest="map1", 
                  default = '30',
                  help="The Mapping length you want!!")
parser.add_option("-M", "--map2", action="store", 
                  dest="map2", 
                  default = '30',
                  help="The 2 Dataset Mapping length you want!!")
parser.add_option("-r", "--readL1", action="store", 
                  dest="readL1", 
                  default = '60',
                  help="The Reads length you want!!")
parser.add_option("-R", "--readL2", action="store", 
                  dest="readL2", 
                  default = '60',
                  help="The 2 Dataset Reads length you want!!")

parser.add_option("-d", "--K_MER_FILTER", action="store", 
                  dest="k_mer_freq", 
                  default = '0',
                  help='''remove low-frequency K-mers with frequency no larger than [default 0] ''')

parser.add_option("-D", "--Edge_TRIM", action="store", 
                  dest="edge_trim", 
                  default = '1',
                  help='''remove edges with coverage no larger than [default 1] ''')

parser.add_option("-x", "--Similar_merge", action="store", 
                  dest="similar", 
                  default = '1',
                  help='''strength of merging similar sequences during contiging [default 1, min 0, max 3]''')

parser.add_option("-f", "--FIX_inner_gap", action="store_true", 
                  dest="fix", default= False,
                  help='''intra-scaffold gap closure''')

parser.add_option("-g", "--Allow_estimate", action="store", 
                  dest="deviation", 
                  default = '50',
                  help='''allowed length difference between estimated and filled gap[ default, 50  ]''')



parser.add_option("-L", "--Allow_contig", action="store", 
                  dest="allow_contig", 
                  default = '100',
                  help='''minimum contigs length used for scaffolding''')

(options, args) = parser.parse_args() 
'''---------------The paramater add new ----------------------'''
k_mer_freq = options.k_mer_freq

edge_trim = options.edge_trim

similar = options.similar

fix = options.fix

deviation = options.deviation

allow_contig = options.allow_contig

'''-------------------------------------------------------------------'''

s1      = options.s1
s2      = options.s2

cpu     = options.cpu
read1   = options.read1
read2   = options.read2
read3   = options.read3
read4   = options.read4
thread  = options.thread
inszie1 = options.inst1
inszie2 = options.inst2
map_length1 = options.map1
map_length2 = options.map2
read_length1 = options.readL1
read_length2 = options.readL2
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
if __name__=='__main__':
	for i in need_r:
		CONFIG = open('./K%s/config%s.ini'%( i,i ),'w')
		CONFIG.write('''max_rd_len=100\n[LIB]\navg_ins=%s\nreverse_seq=%s\nasm_flags=3\nrank=1\npair_num_cutoff=5\nrd_len_cutoff=%s\nmap_len=%s\nq1=%s\nq2=%s\n[LIB]\navg_ins=%s\nreverse_seq=%s\nasm_flags=3\nrank=2\npair_num_cutoff=5\nrd_len_cutoff=%s\nmap_len=%s\nq1=%s\nq2=%s\n'''%( inszie1,s1,read_length1,map_length1,read1,read2,inszie2,s2,read_length2,map_length2,read3,read4 ))
	def run(x):
		x=str(x)
		if fix:
			command = 'soapdenovo63mer all -s ./K'+x+'/config'+x+'.ini -K'+x+\
				    ' -R -o ./K'+x+'/soap_k_'+x+' -p '+ str( cpu )+ \
				    ' -d ' +k_mer_freq +' -D '+edge_trim+' -M '+similar+' -F '+' -G '+deviation+' -L '+allow_contig   \
				    +' >run'+x+'.log '
		else:
			command = 'soapdenovo63mer all -s ./K'+x+'/config'+x+'.ini -K'+x+\
				    ' -R -o ./K'+x+'/soap_k_'+x+' -p '+ str( cpu )+ \
				    ' -d ' +k_mer_freq +' -D '+edge_trim+' -M '+similar  + ' -L ' + allow_contig   \
				    +' >run'+x+'.log '
		os.system( command)
	print( need_r )	
	pool = multiprocessing.Pool(thread)
	pool.map( run, need_r)

	END = open( 'Assembly_status.output','w' )
	for e_f in glob.glob( '*.log' ):
		try:
			kmer = re.search( 'run(\d+)\.log',e_f ).group(1)
		except:
			continue
		
		RAW = open( e_f,'rU' )
		data = re.search( '\n(\d+\s+scaffolds from .+\n.+\n.+\n.+\n)' ,RAW.read() )
		if data:
			END.write( '%s\n>>>>>\n%s-----------------------------\n'%( kmer,data.group(1) )  )