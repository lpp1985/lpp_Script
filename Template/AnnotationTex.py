#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/11/12
"""


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
	template = env.get_template('Annotation.tex')
	END = open( InputPath+"Annotation.tex" ,'w' )

	result_dir = InputPath+"04.OtherDatabase/"
	all_cog_tex = []
	all_kegg_tex = []
	for a,b,c in os.walk( result_dir+"Detail/" ):
		name = a.split("/")[-2].split("_")[-1]
		if a.endswith("eggNOG"):
			if os.path.exists(a+"/stats.pdf"):
				result = a+"/stats.tex"
				commandline = """Graph2tex.py  -i %s  -o %s -c %sCOG分布统计柱状视图  """%(    
				a+"/stats.pdf",result,name
			) 
				
				os.system( commandline  )
				all_cog_tex.append(result)
		if a.endswith("KEGG"):
			
			result = a+"/stats.tex"
			if os.path.exists(a+"/stats.pdf"):
				commandline = """Graph2tex.py  -i %s  -o %s -c %sKEGG通路统计柱状视图  """%(    
					        a+"/stats.pdf",result,name
					    ) 
			
				all_kegg_tex.append(result)			
				os.system( commandline )

	
	all_kegg_tex = sorted( all_kegg_tex,key = lambda x: os.path.dirname( x ).split("/")[-2]  )
	all_cog_tex = sorted( all_cog_tex,key = lambda x: os.path.dirname( x ).split("/")[-2]  )
	table_dir = result_dir+"Table/Database"
	commandline =  """ txt2latex.py -i %s/stat.tsv -o %s/stat.tex   -c  "注释结果统计表" """%(  table_dir,table_dir  )  

	os.system( commandline   )    
	table = "%s/stat.tex"%( table_dir  )
	commandline =  """Graph2tex.py -i %s/stat.pdf -o %s/stat_graph.tex   -c  "注释结果数据库来源分布图" """%(  table_dir,table_dir  )  
	graph = "%s/stat_graph.tex"%( table_dir  )
	os.system( commandline)  
	
	END.write(
	    template.render(
	        {
	            "root_path":template_root,
	            "COGGRAPH":all_cog_tex,
	            "kegggraph":all_kegg_tex,
	            "annotationtable":table,
	            "annotationgraph":graph,
	        }
	        ).encode('utf-8')
	)	
	END.close()	


