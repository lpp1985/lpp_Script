#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/1/23
"""
from optparse import OptionParser
import sys,os
from GO_obo_parse import *
usage = '''usage: python2.7 %prog -i input_path -t [The type you want]'''
parser = OptionParser(usage =usage ) 
parser.add_option("-i", "--INPUT", action="store", 
                  dest="input",
                  help="Input File")

parser.add_option("-o", "--end", action="store", 
                  dest="output", 
                  help="OUTPUT Data")
if __name__=="__main__":
	(options, args) = parser.parse_args() 
	RAW = open(options.input,'rU')
	END = open(options.output,'w')
	RAW.next()
	END.write("Component\tFunction\tID\tGeneNumber\n")
	data_hash = Ddict()
	for line in RAW:
		line_l = line.split("\t")
		if line.startswith("GO:"):
			function =line_l[1]
			goid = line_l[0]
			component =  GO_COMPONENT.select(GO_COMPONENT.q.Go==line_l[0])[0].Compent
			component = component.capitalize().replace("_"," ")
		else:
			gene = line_l[3]
			data_hash[component][function+'\t'+goid][gene] = ""
	
	for compe,func_hash in data_hash.items():
		for func ,gene_hash in func_hash.items():
			
			END.write(compe+'\t'+func+'\t%s\n'%(len(gene_hash)))

		
