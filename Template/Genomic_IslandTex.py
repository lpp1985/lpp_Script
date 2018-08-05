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
	template = env.get_template('Genomic_Island.tex')


	result_dir = InputPath+"08.Genomic_Island"
	all_IS = []

	for a,b,c in os.walk( result_dir):
		name = a.split("/")[-1].split("_")[-1]
		a = a+'/'
		for f in c:
			if f.endswith("zb195_Genome2_GI.stat"):

				result = a+"/stats.tex"
				data = pd.read_table(a+f)
				if len(data)==0:
					continue

				os.remove(  a+f+".matrix" )
				commandline = """txt2latex.py  -i %s  -o %s -c "%s  基因岛分类统计表"  """%(    
				    a+f,result,name
				) 

				os.system( commandline  )
				all_IS.append(result)


	all_IS = sorted( all_IS,key = lambda x: os.path.dirname( x ).split("/")[-1]  )


	END = open( InputPath+"Genomic_Island.tex" ,'w' )
	END.write(
	    template.render(
	        {
	            "GI":all_IS,

	        }
	        ).encode('utf-8')
	)	
	END.close()	


