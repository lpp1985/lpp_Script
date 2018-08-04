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

parser.add_option("-r", "--REFERENCE", action="store",
		              dest="reference_all",

		              help="The reference appendix name to blast")

parser.add_option("-i", "--ISLAND", action="store",
                              dest="island",

                              help="The island file")

parser.add_option("-o", "--OUTPUT", action="store",
                              dest="output",

                              help="The output Result")
parser.add_option("-a", "--All_Gene", action="store",
                              dest="all_gene",

                              help="All the gene file ( fasta format  )")
(options, args) = parser.parse_args()
all_gene = options.all_gene



OUTPUT = open( options.output,'w'   )

ALL_GENE = fasta_check( open( all_gene ,'rU'   ) )

island = options.island

ISLAND = fasta_check( open( options.island ,'rU'   ) )

all_reference = glob.glob( '*.'+options.reference_all  )
total_result = []
def solve_coordinate(  line  ):
		
	
	line_l = line.strip().split('\t')
	all_perc = float( line_l[-5] ) 
	
	if all_perc <=90 :
		return ''
	if  line_l[-1] not in number_hash:
		number_hash[ line_l[-1] ] =''
	else:
		
		no_unique[ line_l[-1] ]	 = ''
def solve_island_coordinate(  line  ):

	global no_range,all_island_location_hash
	
	line_l = line.strip().split('\t')
	

	all_perc = float( line_l[-5] ) 
	
	if all_perc <=90 :
		return ''
	location_list = map( int, [ line_l[0],line_l[1] ] )		
	all_island_location_hash+= [location_list[0] ,location_list[1]]
	
	no_range+=xrange(  location_list[0] ,location_list[1] )

def solve_island_near_coordinate(  line  ):

	global have
	line_l = line.strip().split('\t')
	all_perc = float( line_l[-5] ) 
	
	if all_perc <=90 or line_l[-1] not in unique :
		return ''
	
	location_start = int(line_l[0] )
	location_stop = int( line_l[1]  )
	gene_location = (location_start +location_stop ) /2
	
	if gene_location  in no_range:
		return ''

	for each_location in sorted( all_island_location_hash ):
		if abs( each_location - gene_location ) <= 20000:

			have[ line_l[-1]  ] = ''
for each_reference in all_reference:
	
	
		
	no_range = []
	have = {}	
	
	all_island_location_hash = []
	number_hash = {}
	no_unique = {}
	unique = {}
	name = os.path.split( each_reference  )[-1].split('.')[0]
	os.system( 'nucmer --maxmatch %s %s 2>/dev/null'%( each_reference , all_gene  ) )
	#os.system( 'delta-filter  -rqg out.delta  >out2.delta' )
	os.system( 'show-coords  -cdTHr out.delta  >%s.cache'%(  name  ) )	

	RAW = open( '%s.cache'%( name  ) ,'rU' )
	map( solve_coordinate ,RAW )
	for each_data in number_hash:
		if each_data not in no_unique:
			unique[ each_data ] = ''
	
	os.system( 'nucmer --maxmatch %s %s 2>/dev/null'%( each_reference , island  ) )
	#os.system( 'delta-filter  -rqg out.delta  >out2.delta' )
	os.system( 'show-coords  -cdTHr out.delta  >%s.cache_island'%(  name  ) )
	RAW = open( '%s.cache_island'%( name  ) ,'rU' )
	
	map( solve_island_coordinate ,RAW )
	
	RAW = open( '%s.cache'%( name  ) ,'rU' )
	
	map(  solve_island_near_coordinate ,RAW    )

	total_result.append ( set( have  )  )



def get_detail(SEQ):
	end_hash = {}
	for t,s in SEQ:
		end_hash[  t[1:-1].split()[0]  ] = re.sub( '\s+','',s )
	return end_hash

all_gene_seq_hash = get_detail( ALL_GENE  ) 



	
	
	
if __name__=='__main__':
	all_have_gene = reduce( lambda x,y: x&y , total_result  )
	for each_f in all_have_gene:
		print( each_f )
