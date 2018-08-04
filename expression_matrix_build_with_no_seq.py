#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/6/2
from lpp import *
#usage python2.7 expression_matrix_build.py     [ count.py's end__1_append_file_name  ]  [   matrix_end   ] [ FASTA_FILE ]
import glob,re,sys
all_f = glob.glob('*.%s'%(sys.argv[ 1 ] )  )
list_cand = sorted(all_f)
ALL_FASTA = fasta_check(   open( sys.argv[2], 'rU'  )    )
all_seq = {}
for t,s in ALL_FASTA:
	all_seq[ t[1:-1] ] = t+s
for files_1 in list_cand:
	for files_2 in list_cand[  list_cand.index( files_1 )  :]:
		if files_1!=files_2:
			name_1 = files_1.split( '.' )[0]; name_2=files_2.split( '.' )[0]
			END = open( name_1+'_'+ name_2+'.matrix','w'  )
			STATIC = open(  name_1+'_'+ name_2+'.check' ,'w')
			UNIQUE = open(  name_1+'_'+ name_2+'.unique' ,'w')
			NONE = open(  name_1+'_'+ name_2+'.none' ,'w')
			data1 = re.sub( '\t(\S+)',lambda x: '\t'+str(int( float( x.group(1) ) )  )   ,open( files_1  ).read()  )
			data2 = re.sub( '\t(\S+)',lambda x: '\t'+str(int( float( x.group(1) ))  )   ,open( files_2  ).read()  )
			d1_hash = dict(re.findall(  '(\S+)\t(\S+)',data1 ))
			d2_hash = dict(re.findall(  '(\S+)\t(\S+)',data2 ))
			all_have = {}
			END.write(  'gene\t%s\t%s\n'%( name_1,name_2 )  )
			STATIC.write(   'gene\t%s\t%s\n'%( name_1,name_2 )  )
			unique1 = 0
			unique2 = 0
			for key1 in d1_hash:
				if key1 in d2_hash:
					end = d2_hash[key1]
				else:  
					end = '0'
				if end=='0':

					STATIC.write( '%s\t%s\t%s\n'%( key1, d1_hash[key1] ,end  ) )
					unique1+=1
				else:
					
					END.write( '%s\t%s\t%s\n'%( key1, d1_hash[key1] ,end  )  )
				all_have[key1] = ''
			for key2 in d2_hash:
				if key2 not in d1_hash:

					unique2+=1
					#END.write( '%s\t0\t%s\n'%( key2, d2_hash[key2]  )  )
					STATIC.write( '%s\t0\t%s\n'%( key2, d2_hash[key2]  )  )
				all_have[key2] = ''
			UNIQUE.write(  '%s\t%s\n'%( name_1,name_2 )  )
			UNIQUE.write( '%s\t%s\n'%(unique1, unique2) )
			ALL_HAVE = open( name_1+'_'+ name_2+'.validate','w'   )
			ALL_HAVE.write( '\n'.join(  all_have  )+'\n' )
			FASTA = open(  name_1+'_'+ name_2+'.fasta', 'w' )
			for each_t in all_have:
				FASTA.write( all_seq[each_t ]  )