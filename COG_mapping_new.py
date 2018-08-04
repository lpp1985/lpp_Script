#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/19
from lpp import *
from optparse import OptionParser 
from parse_eggNog import *
usage = '''usage: python2.7 %prog -i input_path -t [The type you want]'''
parser = OptionParser(usage =usage ) 
parser.add_option("-i", "--INPUT", action="store", 
                  dest="input",
                  default = './', 
                  help="Input File")

parser.add_option("-c", "--COG", action="store", 
                  dest="cog",
                  default = 'COG', 
                  help="COG,NOG or KOG")

parser.add_option("-o", "--end", action="store", 
                  dest="output", 
                  help="OUTPUT Data")
(options, args) = parser.parse_args() 

BLAST = open(options.input,'rU')
BLAST.next()
END = open(options.output,'w')
END.write("Gene\tSubject\tE_value\tBitscore\tCOG\tAnnotation\tCat\tCategory Annotation\n")
for line in open(options.input,'rU'):
	line_l = line.strip().split("\t")
	subj= line_l[1].split()[0]
	score = line_l[-1]
	e_value = line_l[-2]
	Bit_Score = line_l[-1]
	Identity = line_l[2]
	StartQuery = line_l[7]
	Mismatch = line_l[5]
	EndSubj = line_l[9]
	StartSubj = line_l[8]	
	EndQuery = line_l[8]
	query = line_l[0].split()[0]
	line_l[0] = query
	gene_nog = NOG_GENE.select(AND(NOG_GENE.q.Gene==subj, NOG_GENE.q.NOG.startswith(options.cog))   )
	unique = {}
	for each_gene_nog in gene_nog:
		
		description = NOG_des.select(NOG_des.q.Name==each_gene_nog.NOG)[0].Description
		nog_cat = [NOG_CAT.select( NOG_CAT.q.NOG==each_gene_nog.NOG  )[0]]
		for each_cat in nog_cat:
			cat_anno = CAT_DES.select(CAT_DES.q.Abb==each_cat.Cat  )
			for each_anno in cat_anno:
				END.write("\t".join(line_l)+"\t"+each_gene_nog.NOG+'\t'+description+'\t'+each_cat.Cat+'\t'+each_anno.Description+'\n')
	
	
	
