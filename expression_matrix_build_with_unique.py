#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/6/3
from lpp import *
#usage python2.7 expression_matrix_build.py     [ count.py's end__1_append_file_name  ]  [   matrix_end   ]
import glob,re,sys
all_f = sys.argv[1:]
list_cand = sorted(all_f)
for files_1 in list_cand:
	for files_2 in list_cand[  list_cand.index( files_1 )  :]:
		if files_1!=files_2:
			name_1 = files_1.split( '.' )[0]; name_2=files_2.split( '.' )[0]
			END = open( name_1+'_'+ name_2+'.matrix','w'  )
			data1 = re.sub( '\t(\S+)',lambda x: '\t'+str(int( float( x.group(1) ) )  )   ,open( files_1  ).read()  )
			data2 = re.sub( '\t(\S+)',lambda x: '\t'+str(int( float( x.group(1) ))  )   ,open( files_2  ).read()  )
			d1_hash = dict(re.findall(  '(\S+)\t(\S+)',data1 ))
			d2_hash = dict(re.findall(  '(\S+)\t(\S+)',data2 ))
			
			END.write(  'gene\t%s\t%s\n'%( name_1,name_2 )  )
			U1_unique = open(  name_1+'UniquE__'+name_2+'.unique'  ,'w' )
			for key1 in d1_hash:
				if key1 in d2_hash:
					end = d2_hash[key1]
					END.write( '%s\t%s\t%s\n'%( key1, d1_hash[key1] ,end  )  )
					
				else:
					U1_unique.write( '%s\t%s\t%s\n'%( key1, d1_hash[key1] ,'0'  )  )
					
					END.write( '%s\t%s\t%s\n'%( key1, d1_hash[key1] ,'0'  )  )
					
			U2_unique = open(  name_1+'__UniquE'+name_2+'.unique'  ,'w' )	
			for key2 in d2_hash:
				if key2 not in d1_hash:
					U2_unique.write( '%s\t0\t%s\n'%( key2, d2_hash[key2]  )  )
					END.write( '%s\t0\t%s\n'%( key2, d2_hash[key2]  )  )