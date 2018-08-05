#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/6/3
from lpp import *
#usage python2.7 expression_matrix_build.py     [ count.py's end__1_append_file_name  ]  [   matrix_end   ]
import glob,re,sys
all_f = sys.argv[1:]
total_data = Ddict()
output_name = []
for each_f in all_f:
	
	name = os.path.split(os.path.abspath(each_f))[-1].replace(".count","")
	output_name.append( name )
	total_data[ name ] = File_dict( open( each_f,'rU' )  ).read(1,2)
all_title = {}
for each_key in total_data:
	all_title.update( total_data[ each_key ] )
path = os.path.abspath( os.path.split(sys.argv[1])[0]  )+'/'	
END = open( path+ '__'.join( output_name )+'.matrix','w' )

END.write( 'gene' )
for each_name in output_name:
	END.write( '\t'+each_name )
END.write( '\n' )
for each_gene in all_title:
	END.write(each_gene)
	for each_sample in output_name:
		if each_gene not in total_data[ each_sample ]:
			express = '\t0'
		else:
			express = '\t%s'%(  total_data[ each_sample ][ each_gene ] )
		END.write( express )
	END.write( '\n' )
