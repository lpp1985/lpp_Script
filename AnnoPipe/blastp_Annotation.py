#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/12/6
"""
from optparse import OptionParser
import os
import pandas as pd
from ConfigParser import ConfigParser


if __name__=="__main__":
	usage = '''usage: python2.7 %prog'''
	parser = OptionParser(usage =usage ) 
	parser.add_option("-i", "--INPUT", action="store", 
                      dest="input", 
                      help="input file")

	parser.add_option("-o", "--end", action="store", 
                      dest="output", 
                      help="output")
	parser.add_option("-d", "--Database", action="store", 
	                  dest="database", 
	                  help="database")	

	parser.add_option("-e", "--evalue", action="store", 
                      dest="evalue", 
                      help="evalue cutoff")
	(options, args) = parser.parse_args() 


	Database = options.database
	temp_name = str(os.getpid())
	base_path = os.path.split(os.path.abspath(options.output))[0]+'/'
	if not os.path.exists(base_path):
		os.makedirs(base_path)
	tmp = base_path+temp_name
	os.system(""" blastp -db %s -query %s  -num_threads 64 -max_target_seqs 1 -evalue %s  -outfmt  5    >> %s"""% ( Database,options.input,options.evalue,tmp+'.xml'  )   )
	os.system( "blast_parse.py %s"%(tmp+'.xml') )
	os.remove(tmp+'.xml')
	if os.path.getsize(tmp+".Bparse"):
		os.system(  "cut_best1.py  -i %s.Bparse -o %s.top1  -f "%( tmp,tmp  )  )
		os.remove( "%s.Bparse"%(tmp)  )
		
		dataframe = pd.read_table(  tmp+'.top1'  )
		dataframe = dataframe.drop("Nt_num",1  )
		dataframe = dataframe.drop("Nt_num.1",1  )
		dataframe = dataframe.drop("Nt_score",1  )
		dataframe = dataframe.drop("Nt_query-frame",1  )  
		dataframe = dataframe.drop("Nt_positive",1  )
		dataframe = dataframe.drop("Nt_align-len",1  ) 
		dataframe["Nt_Hit"] = dataframe["Nt_id"].astype("string")+" " + dataframe["Nt_def"].astype("string")
		dataframe = dataframe.drop(["Nt_id","Nt_def","Nt_accession"],1  ) 
		dataframe.rename( columns={"Nt_len":"Nt_Hit_len"},inplace=1     )
		dataframe["Nt_Hit_Coverage"] = '('+(  abs(  dataframe["Nt_hit-from"]-dataframe["Nt_hit-to"]  )+1  ).astype("string")+'/'+dataframe["Nt_Hit_len"].astype("string")  +') '+ ( (  abs(  dataframe["Nt_hit-from"]-dataframe["Nt_hit-to"]  )+1  )/dataframe["Nt_Hit_len"] ).astype("string")
		
		dataframe["Nt_Query_Coverage"] = '('+(abs(  dataframe["Nt_query-from"]-dataframe["Nt_query-to"]  )+1).astype("string")+'/'+dataframe["Nt_query-len"].astype("string")  +') '+ (100*(  abs(  dataframe["Nt_query-from"]-dataframe["Nt_query-to"]  )+1  )/dataframe["Nt_query-len"]).astype("string")
		
		dataframe["Nt_identity"]  = 100*dataframe["Nt_identity"]/dataframe["Nt_query-len"]
		dataframe2 = pd.DataFrame(dataframe,columns=['Name',u'Nt_Hit',"Nt_identity", 'Nt_query-len',u'Nt_Query_Coverage', u'Nt_query-from',u'Nt_query-to', u'Nt_Hit_len','Nt_Hit_Coverage',u'Nt_hit-frame',  u'Nt_hit-from',u'Nt_hit-to', u'Nt_evalue' , u'Nt_bit-score' ])
		dataframe2.to_csv(options.output,sep="\t",index=False)
		os.remove( "%s.top1"%(tmp)  )


	else:
		os.remove(tmp+".Bparse")
		END = open(options.output,'w')