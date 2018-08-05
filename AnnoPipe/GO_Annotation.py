#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/12/6
"""
from Dependcy import *
from optparse import OptionParser

if __name__=="__main__":
	usage = '''usage: python2.7 %prog'''
	parser = OptionParser(usage =usage ) 
	parser.add_option("-i", "--INPUT", action="store", 
	                  dest="input", 
	                  help="input file")

	parser.add_option("-o", "--end", action="store", 
	                  dest="output_prefix", 
	                  help="output_prefix")

	parser.add_option("-e", "--evalue", action="store", 
	                  dest="evalue", 
	                  help="evalue cutoff")





	(options, args) = parser.parse_args() 
	FASTA = fasta_check(open(options.input,'rU'))
	sequence = FASTA.next()[-1]
	blast_type = Nul_or_Protein(sequence)
	output_prefix = os.path.abspath(  options.output_prefix )
	out_put_path = os.path.split(output_prefix)[0]+'/'
	

	if not os.path.exists( out_put_path ):
		os.makedirs( out_put_path )
	if os.path.exists(out_put_path+'/stats.pdf'):
		sys.exit()
	README = open(out_put_path+"Readme.txt",'w')
	README.write(
"""
将所有的基因序列比对到swissprot数据库，并使用blast2go进行GO Mapping。所有的GO根据GO的有向无环图向上回溯，直到第三层。而后统计每一个第三层
结果说明如下：
*.GO-mapping.detail\t基因映射到GO的列表，可以提交到WEGO等网站自动生成可视化结果。EXCEL打开
*.Genome1.GO-mapping.list\t根据有向无环图自动回溯的GO过程，包含每一个第三级GO下所包含的基因和期GO映射,用Excel打开。
*_GO.stats\t每一个第三级GO所映射到的基因个数，由于GO是有向无环图结果，一个子GO可能具有多个parent GO。所以该部分的基因总数要大于实际的基因总数。Excel打开
*_GO.tsv\t每一个基因的GO详细映射结果，用Excel打开
*_SwissAlignment.tsv\t所有基因与swissprot数据库比对的详细结果，用Excel打开
*.xls\tSwissprot和Gene Ontology分析结果的整合结果，用Excel打开。
stat*\tGO分析的可视化结果。
Draw.R\tGO分析可视化画图脚本，用R运行。

"""	
	
	
	    )
	diamond_result = output_prefix+'_Alignment.tsv'
	error = RunDiamond(options.input,options.evalue, blast_type,"uniprot",diamond_result)
	if error:
		print( colored("%s 's Uniprot process in Diamond of eggnog is error!!","red") )
		print(colored( error,"blue"  ))
		print(  "##############################################"   )

		sys.exit()
	uniprot_anno_frame = pd.read_table(  diamond_result  )	
	
	
	
	gomapping_command = "GO_Mapping.py    -i %s  -o %s"%( diamond_result,output_prefix)
	gomapping_process = subprocess.Popen( gomapping_command.split(),stderr= subprocess.PIPE,stdout=  subprocess.PIPE  )
	stdout,stderr = gomapping_process.communicate()	

	
	getGO_command = "get_GO.py  -i %(result)s.GO-mapping.detail -o   %(result)s_GO.xls"%({ "result":output_prefix  })
	getGO_command_list = getGO_command.split()
	getGO_process = subprocess.Popen( getGO_command_list,stderr= subprocess.PIPE,stdout=  subprocess.PIPE  )
	stdout,stderr = getGO_process.communicate()	

	
	golist_command = "GO_List.py -i %(out)s.GO-mapping.list -o %(out)s_GO.stats"%( {"out":output_prefix}  )
	golist_command_list = golist_command.split()
	golist_process = subprocess.Popen( golist_command_list,stderr= subprocess.PIPE,stdout=  subprocess.PIPE  )
	stdout,stderr = golist_process.communicate()	
	
	
	godraw_command = "GO_Draw.py   -i %s_GO.stats  -o %s -r %s"%(
	    output_prefix,
	    out_put_path+"stats",
	    out_put_path+'Draw.R',
	)
	cogdraw_process = subprocess.Popen(  godraw_command.split(),stderr= subprocess.PIPE,stdout=  subprocess.PIPE  )
	stdout,stderr = cogdraw_process.communicate()	
	os.remove( diamond_result )
	
	






