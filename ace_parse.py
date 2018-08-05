#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/4/13
from lpp import *
from optparse import OptionParser 
usage = '''usage: python2.7 %prog -s sequence_result( 454,or 3730  ) -l gap_location -b blast_result'''
parser = OptionParser(usage =usage )
ok_closure = {}
parser.add_option("-s", "--SEQUENCE", action="store", 
                  dest="seq", 
                  help="Contig_fasta")
parser.add_option("-a", "--ACE", action="store", 
                  dest="ace", 
                  help="ACE  File")
parser.add_option("-n", "--NEED", action="store", 
                  dest="need", 
                  help="The tag you need")
parser.add_option("-o", "--Output", action="store", 
                  dest="output", 
                  help="output file")
(options, args) = parser.parse_args() 
seq_file   = options.seq
SEQ = fasta_check( open( seq_file ,'rU'  )  )

ace  = options.ace
ACE = open( ace ,'rU'  )
ace_list = re.split( '\n(?=CO Contig\d+)',ACE.read() )
output = options.output

need  = re.escape(options.need)
OUTPUT = open( output,'w' )
SEQ = fasta_check( open( seq_file,'rU'  ) )
all_need={}
doubt={}
for each_block in ace_list:
	
	need_list = re.findall( '\nAF ('+need+')',each_block )
	
	if need_list:
		name = re.search('^CO (Contig\d+)',each_block).group(1)
		if len(need_list ) !=1:
			doubt[ name  ] = ''
		
			all_need[ name ] = ''
for t,s in SEQ:
	t_name = re.search( '\.(Contig\d+)\n',t ).group(1)
	#if len(s)>=1000:
	if t_name in all_need:
		OUTPUT.write(t+s)
	#else:
		#if t_name not in doubt:
			#OUTPUT.write(t+s)

		