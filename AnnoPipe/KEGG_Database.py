#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/3/16
"""

#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/3/16
"""
import os,sys

import pandas as pd
from os.path import abspath
sys.path.append( os.path.split(abspath(__file__))[0]+'/../Lib/' )
from lpp import *
from Dependcy import *
usage = "python2.7 %prog [options]"
parser = OptionParser(usage =usage )
parser.add_option("-d", "--Database", action="store",
                  dest="DB_FILE",

                  help="Database File")


parser.add_option("-a", "--ANNO", action="store",
                  dest="ANNO",

                  help="KEGG Ghostz Aligment Result")


parser.add_option("-p", "--PATHWAY", action="store",
                  dest="PATH",

                  help="KEGG Pathway Detail")

if __name__ == '__main__':
	(options, args) = parser.parse_args()


	
	

	kegg_annotation_frame = pd.read_table(options.ANNO)
	

	# kegg_annotation_frame["Function"] = kegg_annotation_frame["KEGG_Hit"].str.split(' ',1,return_type='frame')[1]
	# column_name = list( kegg_annotation_frame.columns[-1:] )
	# column_n2 = list( kegg_annotation_frame.columns[1:-1] )
	# column_name.extend(column_n2)
	# kegg_annotation_frame = pd.DataFrame(kegg_annotation_frame,columns=column_name)
	
	
	
	PATHWAY_DETAIL = open( options.PATH,'rU')
	pathway_detail_frame = pd.read_table( options.PATH )
#	pathway_detail_frame = pd.DataFrame(pathway_detail_frame, columns=pathway_detail_frame.columns[:-1])
	pathway_detail_frame = pathway_detail_frame.drop_duplicates()
	# pathway_detail_frame.rename(  columns={'Gene_ID':'Name'}, inplace=True)
	new_name = [ "KEGG_"+name for name in pathway_detail_frame.columns[1:]]
	name_transHash =  dict( zip(pathway_detail_frame.columns[1:],new_name) )
	pathway_detail_frame.rename(  columns=name_transHash, inplace=True)
	
	result_data_frame = pd.DataFrame.merge( kegg_annotation_frame,pathway_detail_frame,left_on='Name', right_on='Name', how='outer' )
	result_data_frame.to_csv(options.DB_FILE,sep="\t",index =False)
	
