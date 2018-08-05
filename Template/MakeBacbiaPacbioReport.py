#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/11/12
"""


from lpp import *
import os,subprocess
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

	template_root = config_hash["Location"][  "root" ]+"/Template/"




	templeloader = FileSystemLoader(template_root)
	env = Environment(loader = templeloader)
	template = env.get_template('BacteriaPacbio.tex')


	result_dir = InputPath
	


	END = open( InputPath+"Report.tex" ,'w' )
	END.write(
	    template.render(
	        {
	            "RootPath":template_root,

	        }
	        ).encode('utf-8')
	)	
	END.close()	

	os.system( 
	    "QCTex.py -i %s  -o %s  -s %s -t %s -g %s  -n %s" %( InputPath,OutputPath,Sample,Table ,Graph ,Cell  )
	            
	            
	            )
	
	os.system(  "MakeTitleTex.py -o %s -s %s"%( InputPath, Sample  )   )
	os.system(  "Assembly_StatTex.py  -i %s"%( InputPath )  )
	os.system(  "GenePredictionTEX.py -i %s"%(InputPath) )
	os.system(  "Genomic_IslandTex.py -i %s"%(InputPath) )
	os.system(  "IS_Tex.py -i %s"%(InputPath) )
	os.system(  "ProphageTex.py -i %s"%(InputPath) )
	os.system(  "Repeat_Tex.py -i %s"%(InputPath) )
	os.system(  "AnnotationTex.py -i %s"%(InputPath) )
	os.system(  "Crispr_Tex.py -i %s"%(InputPath) )
	os.system(  "xelatex %s"%( END.name  )   )

	os.system(  "cd %s &&bibtex Report"%(   OutputPath   )   )
	
	os.system(  "xelatex %s"%( END.name  )   )
	
