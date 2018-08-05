#!/usr/bin/env python
#coding:utf-8
# Author:  LPP
# Purpose: 
# Created: 2011/11/7
# Company: Chinese Human Genomic Center at Beijing
'''需要装Mummer'''
from lpp import *
from  optparse import OptionParser
usage = '''usage: python2.7 %prog [options] 

to 
'''
parser = OptionParser(usage =usage )

parser.add_option("-s", "--SCAFFOLD", action="store",
		              dest="scaff",

		              help="The scaffold to fill the gap")

parser.add_option("-p", "--PCR", action="store",
                              dest="pcr",

                              help="The PCR Result")

parser.add_option("-o", "--OUTPUT", action="store",
                              dest="output",

                              help="The output Result")
parser.add_option("-l", "--LEFT", action="store",
                              dest="left",

                              help="The left PCR result")
(options, args) = parser.parse_args()
scaff = options.scaff


pcr = options.pcr
OUTPUT = open( options.output,'w'   )

LEFT = open( options.left,'w'   )
PCR = fasta_check( open( options.pcr ,'rU'   ) )
SCAFFOLD = fasta_check( open( options.scaff ,'rU'   ) )
def get_detail(SEQ):
	end_hash = {}
	for t,s in SEQ:
		end_hash[  t[1:-1].split()[0]  ] = re.sub( '\s+','',s )
	return end_hash

pcr_reads_hash = get_detail( PCR  ) 
scaff_reads_hash = get_detail( SCAFFOLD  )
subtition_hash = {}
def solve_coordinate(    ):
	global scaff_reads_hash,subtition_hash
	all_perc =sum(  [ float( x[-5] )    for x in cache     ] )
	
	if all_perc <=40 or len( cache )<2:
		return ''
	direction = dir_hash.keys()[0]
	all_reference_coor = []
	for x in cache:
		all_reference_coor.append( int(x[0])  )
		all_reference_coor.append( int( x[1] )  )
	all_query_coor = []
	
	for x in cache:
			all_query_coor.append( int(x[2])  )
			all_query_coor.append( int(x[3]) )	
	reference_location = [ min( all_reference_coor ),max( all_reference_coor  )+1   ]
	
	query_location = [ min( all_query_coor ),max( all_query_coor  )+1   ]
	
	
	ref_id = ref_hash.keys()[0]
	reference_sub_sequence = scaff_reads_hash[  cache[0][ -2 ]  ][  reference_location[0] : reference_location[1] ]
	length = len(  reference_sub_sequence  )
	name = al_hash.keys()[0]
	tag1 = '<%s>'%( name  )
	tag2 = '</%s>'%( name )
	left_length = length-len(tag1+tag2)
	transfer_seq = tag1+'N'*left_length+tag2
	scaff_reads_hash[  cache[0][ -2 ]  ] = scaff_reads_hash[  cache[0][ -2 ]  ][ :min( all_reference_coor ) ] +transfer_seq +scaff_reads_hash[  cache[0][ -2 ]  ][ max( all_reference_coor ): ]

	#scaff_reads_hash[  cache[0][ -2 ]  ] = re.sub( '(^\S{%s})\S{%s}'%( min( all_reference_coor )-1 ,length ),'\\1'+transfer_seq      ,scaff_reads_hash[  cache[0][ -2 ]  ] )
	
	query_sequence = pcr_reads_hash[  name ]  [ query_location[0]:query_location[1]  ]
	if direction =='-1':
		subtition_hash[ name ] = complement( query_sequence )
	else:
		subtition_hash[ name ] = query_sequence 
	
if __name__=='__main__':
	os.system( 'nucmer --maxmatch %s %s '%( options.scaff , options.pcr  ) )
	os.system( 'delta-filter  -rqg out.delta  >out2.delta' )
	os.system( 'show-coords  -cdTHr out2.delta  >see' )
	RAW = open( 'see','rU'  )
	cache = []
	cache2 = []
	al_hash = {}
	ref_hash = {}
	dir_hash = {}
	line_l = RAW.next().strip().split('\t')
	cache.append( line_l  )
	al_hash [ line_l[-1]  ] = ''
	ref_hash [  line_l[-2]  ] = ''
	dir_hash[ line_l[-3] ] = ''
	for line in RAW:
		line_l = line.strip().split('\t')
		name = line_l[-1]
		ref = line_l[-2]
		direction = line_l[-3]

		if name in al_hash and ref in ref_hash and direction in dir_hash:
			cache.append(   line_l   )

		
		else:	
			
			solve_coordinate()
			cache=[]
			cache.append( line_l  )
			al_hash = {}
			ref_hash = {}
			dir_hash = {}			
			al_hash [ line_l[-1]  ] = ''
			ref_hash [  line_l[-2]  ] = ''
			dir_hash[ line_l[-3] ] = ''			
	else:
		solve_coordinate()
		
	CACHE = open( 'CACHE','w'  )
	for eac_scaff in scaff_reads_hash:
		CACHE.write( '>%s\n%s\n'%( eac_scaff , scaff_reads_hash[ eac_scaff ]  ) )
	CACHE = open( 'CACHE','rU'  )
	cache_data = CACHE.read()
	for each_name in subtition_hash:
		cache_data = re.sub( '<%s>N+</%s>'%( each_name,each_name ) ,subtition_hash[  each_name ] , cache_data  )

	#cache_data = re.sub( '(\S{60})','\\1\n',cache_data )
	OUTPUT.write( cache_data )
	for each_reads in pcr_reads_hash:
		if each_reads not in subtition_hash:
			LEFT.write(  '>'+each_reads+'\n'+pcr_reads_hash[ each_reads] +'\n' )
	#print( len(  subtition_hash[ '269' ] ) )