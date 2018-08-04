#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/10/13
from lpp import *
from  optparse import OptionParser
usage = '''usage: python %prog [options] 

It can automaticly extract sequence from merged.gtf produced by cuffmerge'''
parser = OptionParser(usage =usage ) 
parser.add_option("-i", "--Input", action="store", 
                  dest="input_gtf",
                  type='string',  
                  help="the gtf file produced by cufflinks")

parser.add_option("-g", "--Genome", action="store", 
                  dest="Genome",
                  type='string',  
                  help="the fasta file of Genome")


parser.add_option("-o", "--Output", action="store", 
                  dest="output_path",
                  type='string',  
                  help="the dir to output")
(options, args) = parser.parse_args() 

def get (  seq,data_list     ):
	strand = data_list[0]
	coord_list = data_list[-1]
	seq_output = ''
	for each_coord in coord_list:
		seq_output += seq[ each_coord[0]: each_coord[1]  ]
	if strand == '-':
		
		seq_output = complement(  seq_output  )
	
	seq_output = re.sub( '(\w{60})','\\1\n',seq_output   )
		
	return seq_output


GTF = open( options.input_gtf,'rU'    )

GENOME = fasta_check( open( options.Genome  ,'rU' )   )

output_path = options.output_path

if not os.path.exists(  output_path ):
	os.makedirs(  output_path )
	
'''The Output_coord is a hash to store the coordinate and strand of a tag'''

Output_coord = Ddict()

'''ALL_title stores different ID category '''
All_title = {}


for line in GTF:
	
	line_l = line[:-1].split('\t')
	
	chrome = line_l[0]
	
	coords = [ int(line_l[3])-1  , int( line_l[4] )    ]

	strand = line_l[6]
	
	all_detail = re.findall( '(\S+_id)\s+\"(\S+)\"',line_l[-1]    )
	
	#if chrome not in Output_coord  or :
		
	for each_title,data in all_detail:
		All_title[ each_title ] = ''
		if chrome not in Output_coord or each_title not in Output_coord[ chrome ] or data not in Output_coord[ chrome ] [  each_title ]:
			
			Output_coord[  chrome ] [ each_title ][data] = [   strand,[  coords  ]       ]
			
		else:
			Output_coord[  chrome ] [ each_title ][data][-1].append(  coords  )
			
output_file_hash = {}


for each_title in All_title:
	output_file_hash[ each_title  ] = open( output_path+ each_title.replace('_id','')+'.fasta'  ,  'w' )
	
for t,s in GENOME:
	title = t[:-1].split()[0][1:]

	if title in  Output_coord:
		seq = re.sub('\s+','',s)
		
		for each_category in Output_coord[ title ]:
			
			for each_id in Output_coord[ title ][ each_category  ]:

				out_seq = get(   seq,  Output_coord[ title ][ each_category  ][ each_id ]    )
				
				output_file_hash[  each_category  ].write(  '>'+each_id+'\n'+out_seq+'\n'  )
				
			

			
