#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/11/12
"""

from lpp import *
import os
from  jinja2 import FileSystemLoader,Environment
from Dependcy import Config_Parse
from optparse import OptionParser





if __name__ == '__main__':
	usage = '''usage: python2.7 %prog [options] 
'''
	parser = OptionParser(usage =usage )



	parser.add_option("-i", "--InputPath", action="store",
	                  dest="InputPath",

	                  help="Input Path")

	
	(options, args) = parser.parse_args()
	InputPath = os.path.abspath(options.InputPath)+'/'

	config_hash = Config_Parse()

	template_root = config_hash["Location"][  "root" ]+"/Template"
	



	templeloader = FileSystemLoader(template_root)
	env = Environment(loader = templeloader)
	template = env.get_template('Assembly_Stat.tex')
	END = open( InputPath+"/Assembly_Stat.tex" ,'w' )
	
	result_dir = (InputPath+"/01-3.Assembly_END/")
	total = len( glob.glob(  result_dir +"/*.fasta" ) )
	plasmid =len( glob.glob(  result_dir +"/*Plasmid*" ) ) 
	chrom = total-plasmid
	linear = 0
	cir = 0
	for each_file in  glob.glob(  result_dir +"/*.fasta" ):
		RAW = fasta_check(  open(each_file)   )
		for t,s in RAW:
			if "Circle" in t:
				cir +=1
			else:
				linear +=1
	END.write(
	    template.render(
	        {
	            "Cir":cir,
	            "Linear":linear,
	            "Chro":chrom,
	            "Plsmid":plasmid
	        }
	    ).encode('utf-8')
	)	
	END.close()
	
