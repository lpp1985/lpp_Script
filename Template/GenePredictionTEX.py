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


	(options, args) = parser.parse_args()
	InputPath = os.path.abspath(options.InputPath)+'/'

	config_hash = Config_Parse()

	template_root = config_hash["Location"][  "root" ]+"/Template"




	templeloader = FileSystemLoader(template_root)
	env = Environment(loader = templeloader)
	template = env.get_template('GenePrediction.tex')
	END = open( InputPath+"GenePrediction.tex" ,'w' )

	result_dir = InputPath+"10.CircleGraph/"
	all_graph = []
	all_png = glob.glob(result_dir+"/*.png")
	for each_f in all_png:
		if ".thu.png" not in each_f and  each_f+".thu.png" not in all_png:
			os.system(  "thu.py %s"%(each_f)  )
			
	for each_f in glob.glob(result_dir+'/*.png.thu.png'):
		name = os.path.basename( each_f).split(".")[0].split("_")[-1]
		tex = each_f.replace("png.thu.png","tex")
		commandline = """Graph2tex.py  -i %s  -o %s -c %s基因组视图  """%(    
		    each_f,	tex,name
		) 

		os.system( commandline )

		all_graph.append(  tex )
	all_graph = sorted( all_graph,key = lambda x: os.path.basename( x ).split(".")[0].split("_")[-1]  )
	total_dir = InputPath+"09.AllResult/"
	os.system(  "cd %s && N50-new.py %s"%(  total_dir, "Total.ffn" )    )	
	os.system(  "cd %s && lengthN50.R %s"%(  total_dir, "Total.scope  result.pdf result.tiff" )    )	
	
	commandline = """Graph2tex.py  -i %s/result.pdf  -o %s/result.tex -c 基因长度分布统计图 """%(    
		    total_dir,total_dir
		) 	
	os.system( commandline )
	
	
	lengh_graph = "%s/result.tex"%(  total_dir )                                  
	anno_path = InputPath+"03.Annotation/"
	commandline =  """ txt2latex.py -i %s/stats.tsv -o %s/stats.tex   -c  "序列注释结果统计表" """%(  anno_path,anno_path  )  

	os.system( """ txt2latex.py -i %s/stats.tsv -o %s/stats.tex   -c  "序列注释结果统计表" """%(  anno_path,anno_path  )    )    
	table = "%s/stats.tex"%( anno_path  )
	END.write(
	    template.render(
	        {
	            "CircularGraph":all_graph,
	            "status_table":table,
	            "length_png":lengh_graph,
	         
	        }
	    ).encode('utf-8')
	)	
	END.close()	
	
	
	