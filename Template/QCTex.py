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
	parser.add_option("-o", "--Output", action="store",
	                  dest="OutputPath",

	                  help="Output Path")	
	parser.add_option("-s", "--Sample", action="store",
	                  dest="Sample",

	                  help="Sample Name")		
	
	parser.add_option("-t", "--Table", action="store",
	                  dest="Table",

	                  help="Qc Tex Table")		
	parser.add_option("-g", "--Graph", action="store",
	                  dest="Graph",

	                  help="Qc Graph")			
	
	parser.add_option("-n", "--Number", action="store",
	                  dest="Number",
	                  type=int,
	                  default=1,
	                  help="Cell Number")	
	(options, args) = parser.parse_args()
	InputPath = os.path.abspath(options.InputPath)+'/'
	OutputPath = os.path.abspath(options.OutputPath)+'/'
	Sample = options.Sample
	Cell = options.Number
	Graph = os.path.abspath(options.Graph)
	Table = os.path.abspath(options.Table)
	config_hash = Config_Parse()

	template_root = config_hash["Location"][  "root" ]+"/Template"
	



	templeloader = FileSystemLoader(template_root)
	env = Environment(loader = templeloader)
	template = env.get_template('QC.tex')
	END = open( InputPath+"/QC.tex" ,'w' )
	END.write(
	    template.render(
	        {
	            "Sample":Sample,
	            "Cell":Cell,
	            "Graph":Graph,
	            "Table":Table
	        }
	    ).encode('utf-8')
	)	
	END.close()
	