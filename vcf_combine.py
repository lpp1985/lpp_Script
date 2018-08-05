#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/9/28
from lpp import *
from  optparse import OptionParser
'''a Multiple-division hash to store the data from different VCF'''
all_data_Ddict = Ddict()
'''store all the sample name'''

all_sample = {}

'''store each sample in each SNP pos status'''
sample_status = Ddict()
usage = '''usage: python2.7 %prog -i [Input_PATH] -o [Prefix] -a [ Appendix  ]'''
parser = OptionParser(usage =usage ) 
parser.add_option("-i", "--INPUT", action="store", 
	dest="input_path", 
	help="Input_PATH")
parser.add_option("-o", "--prefix", action="store", 
	dest="prefix", 
	default = 'ALL_vcf',
	help="The PREFIX you want!!!!")
parser.add_option("-a", "--append", action="store", 
	dest="append", 
    default = 'vcf', 
	help="The Thread number  you want !!!!")
(options, args) = parser.parse_args() 

'''The path you want'''
input_path =  options.input_path
'''The Output you like'''
OUTPUT = open(  options.prefix,'w'   )
'''The input apppendix to store the vcf data. and make it trans to regex expression'''
append =  re.escape(  options.append  )
'''traverse the input path and read all the vcf file include : Reference, Pos, REF, ALT,Depth, Varation/RAW,Is it Total identical to all Sample(Yes:*;NO:'-')'''
'''A function to get the data in the vcf file,Input is the file_handle of VCF file ,attention please that the file name must start with sample and follow with.+append'''
def Get__data(  FILE_HANDLE  ):
	isinstance( FILE_HANDLE,file )
	'''Get the sample name ,is the first area of character of file name'''

	Sample_name = os.path.split(FILE_HANDLE.name)[-1].split('.')[0]

	all_sample[ Sample_name ] = ''
	for line in FILE_HANDLE:
		'''Jump the header of VCF file'''
		if line.startswith( '#' ):
			continue
		line_l = line.split('\t')
		reference = line_l[0]
		pos = line_l[1]
		ref_char = line_l[3]
		alt_char = line_l[4]
		'''raw* is the reads number to supports raw, and var* is the reads number to support varaation''' 
		[( raw1, raw2, var1, var2 )]=re.findall( 'DP4=(\d+),(\d+),(\d+),(\d+);',line  )
		total_raw = int(raw1)+int(raw2)
		total_var = int( var1 ) + int( var2 )
		compare_char = str(total_var)+'/'+str(  total_raw  )
		
		
		
		all_data_Ddict[ reference+'\t'+pos+'\t'+ref_char  ][ Sample_name  ] = alt_char+'\t'+\
		              compare_char
		
		sample_status[ reference+'\t'+pos+'\t'+ref_char   ][ alt_char ][ Sample_name  ] = ''
		
		

for a,b,c in os.walk(  input_path ):
	
	for each_file in c:
		
		if re.search( '('+append+')$'  ,each_file):

			Get__data( open(  a+'/'+ each_file ) )
'''Compute_consensus status in all sample'''
consensus_hash = {}
for each_data in sample_status:
	if len( sample_status[  each_data   ]  )==1 and len(    sample_status[  each_data   ][   sample_status[  each_data   ].keys()[0]     ]    ) == len( all_sample ):
		consensus_hash[ each_data ] = ''
		
''' write the header in OUTPUT file'''
OUTPUT.write(  '''#POS: The location of Reference 
#RAW: The Base in the Reference
#Sample: The Sample's  SNP or INDEL base
#Sample__status: The Status of Sample mean [   Number of Reads support Variation  ]/[    Number of Reads cmobat Variation    ]
#Identical   Is all the Sample perform identical SNP to Reference.If so "*" else blank\n'''   )
OUTPUT.write( '#REF\tPOS\tRAW'  )
sample_status = sorted(  sample_status )
for sample in all_sample:
	OUTPUT.write( '\t'+sample+'\t'+sample+'__status'  )

OUTPUT.write( '\tIdentical\n' )

for each_pos in all_data_Ddict:
	OUTPUT.write( each_pos )
	for each_sample in all_sample:
		if each_sample in all_data_Ddict[ each_pos ]:
			OUTPUT.write(  '\t'+ all_data_Ddict[ each_pos ][  each_sample ]  )
		else:
			OUTPUT.write(   '\t-\t-'  )
	if each_pos in consensus_hash:
		OUTPUT.write( '\t'+'*'   )
	else:
		OUTPUT.write( '\t'  )
	
	OUTPUT.write('\n')

