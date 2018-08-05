#!/usr/bin/env python
#coding:utf-8
# Author:   --<>
# Purpose: 
# Created: 2013/5/22


from lpp import *
from optparse import OptionParser
import time
usage='''usage: parse clustal out find snp '''
parser = OptionParser(usage =usage )

parser.add_option("-i", "--Clustal", action="store",
                      dest="clustal",
                      type='string',
                      help="Input Clustal out")

parser.add_option("-n", "--Number", action="store",
                      dest="reference_number",
                      type='int',
                      help="Number of Reference")


(options, args) = parser.parse_args()

RAW = open( options.clustal,'rU'  )
reference_number = options.reference_number
all_have = [  ]

all_data = RAW.read()

print( 'read all' )

data_blocks = re.split( '\n\n+', all_data)[1:]
print(  'split all'  )
name_block = data_blocks[0]
align_data = name_block.split('\n')

i=1
detail  = {}
name_tag = re.search('^(\w+\s+)', align_data[0])
name_end = name_tag.end()
print( 'name ok ' )
i=1
reference_name = ''
for each_block in align_data:
	name_tag = re.search('^(\S+\s+)', each_block)
	if name_tag:
		name = name_tag.group(1).strip()
		aln_b = each_block[name_end:]
		if i !=reference_number:
			detail[ name ] = aln_b
		else:
			reference_name = name
			detail[ 'reference' ] = aln_b		
	else:
		aln_b = each_block[name_end:]
		detail[ 'align'  ] = aln_b
	i+=1
def aln_add( data_string  ):
	
	i=1
	global detail,k
	
	align_data = data_string.split('\n')
	for each_block in align_data:
		name_tag = re.search('^(\S+\s+)', each_block)
		if name_tag:
			name = name_tag.group(1).strip()
			aln_b = each_block[name_end:]
			if i !=reference_number:
				detail[ name ] += aln_b
			else:
				detail[ 'reference' ] += aln_b
		else:
			aln_b = each_block[name_end:]
			detail[ 'align'  ] += aln_b
		i+=1

map( aln_add,data_blocks[1:]  )	

data_diff = re.finditer( '([^\*]+)',  detail[ 'align' ]     )
result = {}
all_other = filter( lambda x: x not in [ 'align','reference'   ] ,detail.keys()   )

snp_length = 0
reference = detail[ 'reference' ]
other = detail[  all_other[0]    ]
def format_check( string ):
	string = re.sub( '\W+','',string )
	if not string:
		string ='.'
	return string
for each_diff in data_diff:
	start,end = each_diff.span()
	
	start_coor = start - len(re.findall( '\W',reference[  :start  ]  ))+1
	differ_block = reference[ start :end  ]
	end_coor = len( re.findall( '\w',differ_block  ) )
	if end_coor ==0:

		start-=1
		start_coor-=1
	else:
		if len( re.findall( '\w',other [start:end]  ) )==0:
			start-=1
	
	result[  start_coor  ] = '%s\t%s\t.\t'%(reference_name, start_coor)+format_check(reference[start:end])+'\t'+ format_check (other [start:end] )   +'\t100.00\tPASS\tDP=100\n'

END = open(  'different.list','w'  )
#for key in sorted(  result ):
#	END.write( '%s\t%s'%( key,result[key] )  )
END2 = open( 'all_diff.list','w' )

END.write( '''##fileformat=VCFv4.1
##fileDate=%s
##reference=%s
##INFO=DP,1,Integer,"Total Depth of Coverage"
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'''%( time.strftime('%Y%m%d', time.localtime(time.time())),reference_name  ) )
for key in sorted(  result ):
	END.write( result[key]  )


