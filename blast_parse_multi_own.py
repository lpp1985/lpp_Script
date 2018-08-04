#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/6/21
from lpp import *
import sys,os,shutil
from optparse import OptionParser
import multiprocessing

usage = '''usage: python %prog [options] 

It can build dataset to input gapmap'''

parser = OptionParser(usage =usage ) 
parser.add_option("-p", "--PEP", action="store", 
                  dest="pep",
                  type='string',  
                  help="Protein sequence file (Fasta)")
parser.add_option("-d", "--DNA", action="store", 
                  dest="dna", 
                  help="nucleocid sequence file(Fasta)")

parser.add_option("-t", "--DIR", action="store", 
                  dest="dir", 
                  help="direct blast result( -m7 output )")

parser.add_option("-r", "--REV", action="store", 
                  dest="rev", 
                  help="reversed blast result( -m7 output )")
(options, args) = parser.parse_args() 
for key1 in [ options.pep, options.dna, options.dir,options.rev     ]:
	if not key1:
		print( '''There's some fault in your paramater , please type the paramater --help to get help'''  )
		sys.exit()
def check_path( path ):
	if os.path.exists(path):
		shutil.rmtree(path)
	os.makedirs( path )
for val in (  './Query/blast/','./Query/dna/','./Query/newout/','./Query/pep/','./Query/out/','./Query/outx/','./Query/newoutx/' ,'./Query/blastx'):
		check_path( val)
def parse_blast(  check,file_name  ):

	RAW = open(options.dir)
	RAW.next()
	Best_hash = Ddict()

	#parse_new = file_name+'.parse'
	#PARSE=blast_parse(  open( file_name ) , open(parse_new,'w') )
	#PARSE.parse()
	#RAW = open(parse_new)
	#RAW.next()
	for line in RAW:
		line_l = line[:-1].split('\t')
		e_val = float(  line_l[12]  )
		if check=='dir':
			
			Best_hash[ line_l[2].split()[0] ] [ e_val ] = line[:-1]
		elif check =='rev':
			Best_hash[ line_l[6].split()[0]] [ e_val ] = line[:-1]
	
	for key1 in Best_hash:
		if len( Best_hash[ key1 ] )>1:
			data = Best_hash[ key1 ][     sorted(  Best_hash[ key1 ]  )[0]       ]
		else:
			data = Best_hash[ key1 ][    Best_hash[ key1 ].keys() [0]       ]
		data_l = data.split('\t')
		if check =='dir':
			
			

			BLAST=open('./Query/blast/'+key1+'.blast','w')
			NEWOUT = open('./Query/newout/'+key1+'.out','w')
			OUT = open('./Query/out/'+key1+'.out','w')
		else:
			BLAST=open('./Query/blastx/'+key1+'.blast','w')
			NEWOUT = open('./Query/newoutx/'+key1+'.out','w')
			OUT = open('./Query/outx/'+key1+'.out','w')
		query = data_l[2]
		query_len = data_l[3]
		sub_len = data_l[20]
		query_align = data_l[ 10 ]
		sub_align = data_l[ 10 ]
		annotation = data_l[ 6 ]
		identity = data_l[ 20 ]+'/'+ query_len
		identiry_rates ='%.4f'%(float( int( data_l[ 20 ] ))/ int(query_len) )
		location = '\t'.join( data_l[13:17] )
		identity = data_l[ 20 ]+'/'+ query_len
		location = '\t'.join( data_l[13:17] )
		score = data_l[ 11 ]
		e_vaul = data_l[12]
		OUT.write( 'Query Name	Query Length	Sbjct Length	Query Alignment	Sbjct Alignment	Annotation	Score	E Value	Identity	Identity_Rate	QueryStart	QueryEnd	SubjectStart	SubjectEnd	\n' )
		OUT.write(  '\t'.join( [query,query_len,sub_len,query_align,sub_align,annotation,score,e_vaul,identity,identiry_rates,location ])+'\n' )
		if check=='dir':
			NEWOUT.write( query+'\t'+annotation.split()[0]+'\t'+ annotation.split(':')[0]+'\t'+e_vaul+'\t'+'1'+'\n')
		else:
			NEWOUT.write( annotation+'\t'+query.split()[0]+'\t'+ query.split(':')[0]+'\t'+e_vaul+'\t'+'1'+'\n')
		BLAST.write( '''Q:\t\t\t%s\nS:\t\t\t%s\nA:\t\t\t%s\n'''%( data_l[-3],data_l[-2],data_l[-1]  )  )
def seq_split( input_file):
	if input_file== options.pep:
		append = '.pep'; path = './Query/pep/'
	else:
		append = '.dna'; path = './Query/dna/'
	RAW= fasta_check(  open(input_file,'rU')  )
	for t,s in RAW:
		END = open( path+t[1:-1].split()[0]+append, 'w' )
		END.write(t+s)
for ele in [ options.pep, options.dna  ]:
	seq_split(  ele )
P_blast  = multiprocessing.Process( target = parse_blast, args = ( 'rev', options.rev )  )
P_blastx = multiprocessing.Process( target = parse_blast, args = ( 'dir', options.dir )  )
P_blast.start()
P_blastx.start()
