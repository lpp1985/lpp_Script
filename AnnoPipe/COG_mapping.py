#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/19
from lpp import *
from optparse import OptionParser 
from parse_eggNog import *
import tempfile
import pandas as pd
usage = '''usage: python2.7 %prog -i input_path -t [The type you want]'''
parser = OptionParser(usage =usage ) 
parser.add_option("-i", "--INPUT", action="store", 
                  dest="input",
                  default = './', 
                  help="Input File")


parser.add_option("-o", "--end", action="store", 
                  dest="output", 
                  help="OUTPUT Data")
(options, args) = parser.parse_args() 
if __name__ == '__main__':
	all_eggnog_frame = pd.read_table(options.input)
	BLAST = open(options.input,'rU')
	BLAST.next()
	TMP = open("%s.tmp"%(os.getpid()),'w')

	
	TMP.write("Name\tCOG\tCOG_Annotation\tCOG_FunCat\tCOG_Category Annotation\n")
	for line in open(options.input,'rU'):
		line_l = line.strip().split("\t")
		subj= line_l[1].split()[0]
		score = line_l[-1]
		e_value = line_l[-2]
		
		query = line_l[0].split()[0]
		gene_nog = NOG_GENE.select(NOG_GENE.q.Gene==subj   )
		unique = {}

		if gene_nog.count()>1:
			gene_nog = gene_nog.limit(1)
		for each_gene_nog in gene_nog:
			
			description = NOG_des.select(NOG_des.q.Name==each_gene_nog.NOG)[0].Description
			nog_cat = [NOG_CAT.select( NOG_CAT.q.NOG==each_gene_nog.NOG  )[0]]
			for each_cat in nog_cat:
				cat_anno = CAT_DES.select(CAT_DES.q.Abb==each_cat.Cat  ).limit(1)
				
				for each_anno in cat_anno:
					TMP.write(line_l[0]+"\t"+each_gene_nog.NOG+'\t'+description+'\t'+each_cat.Cat+'\t'+each_anno.Description.strip()+' [%s]\n'%(each_cat.Cat))
	TMP.close()
	
	
	all_cog_frame = pd.read_table(TMP.name)
	all_cog_frame = all_cog_frame.drop_duplicates()
	result_data_frame = pd.DataFrame.merge( all_eggnog_frame,all_cog_frame,left_on='Name', right_on='Name', how='outer' )
	result_data_frame.to_csv(options.output+".xls",sep="\t",index =False)
	result_stat_frame = pd.DataFrame( result_data_frame,columns=("COG_FunCat","COG_Category Annotation")  )
	result_stat_frame = result_stat_frame[ pd.notnull( result_stat_frame["COG_FunCat"]  )   ]

	result_stat_frame.rename(columns={"COG_FunCat":"Gene"},inplace=True)
	result_stat_frame = result_stat_frame.groupby(["COG_Category Annotation"] ).agg([ 'count'])
	
	result_stat_frame.to_csv(TMP.name,sep = "\t")
	STAT = open(options.output+".stats",'w')
	TMP2 = open(TMP.name,'rU')
	TMP2.next()
	title = TMP2.next()
	title.replace("count","Gene count")
	STAT.write("Category\t"+TMP2.next())
	for line in TMP2:

		name = re.search("\[(\w+)\]$",line.split("\t")[0]).group(1)
		STAT.write( name+'\t'+ line)
	os.remove(TMP.name)
	os.remove(TMP2.name)