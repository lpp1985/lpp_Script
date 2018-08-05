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
	
	parser.add_option("-c", "--COG", action="store", 
		              dest="cog",
		              default = 'COG', 
		              help="COG,NOG or KOG")
	
	
	
	(options, args) = parser.parse_args() 
	FASTA = fasta_check(open(options.input,'rU'))
	sequence = FASTA.next()[-1]
	blast_type = Nul_or_Protein(sequence)
	output_prefix = os.path.abspath(  options.output_prefix )
	out_put_path = os.path.split(output_prefix)[0]+'/'
	end_list = glob.glob(out_put_path+'/stats.pdf')
	if end_list:
		sys.exit()
	cog = options.cog
	if not os.path.exists( out_put_path ):
		os.makedirs( out_put_path )
	README = open(out_put_path+'/Readme.txt','w')
	
	README.write(
"""
将所有的序列比对到eggNOG数据库，并映射到COG。该文件夹结果如下：

*_AlignEggNOG.tsv\t与eggNOG数据库的详细比对结果，Excel打开
*_COG.tsv\t每一个基因的COG映射结果,Excel打开
*_COG.xls\t每个基因与eggNOG比对和COG映射结果整合结果,Excel打开
*_COG.stats\t每一个COG功能分类的基因个数统计表,Excel打开



"""	
	
	)
	diamond_result = output_prefix+'_AlignEggNOG.tsv'
	error = RunDiamond(options.input,options.evalue, blast_type,"eggnog",diamond_result)
	if error:
		print( colored("%s 's COG process in Diamond of eggnog is error!!","red") )
		print(colored( error,"blue"  ))
		print(  "##############################################"   )
		
		sys.exit()
	
	cogmapping_command = "COG_mapping.py  -i %s  -c %s -o %s"%( diamond_result,cog,output_prefix+"_COG")
	cogmapping_process = subprocess.Popen( cogmapping_command.split(),stderr= subprocess.PIPE,stdout=  subprocess.PIPE  )
	stdout,stderr = cogmapping_process.communicate()	
	
	
	cogdraw_command = "COG_Draw.py   -i %s.xls  -o %s -r %s"%(
	    output_prefix+"_COG",
	    out_put_path+"stats",
	    out_put_path+'Draw.R',
	)
	cogdraw_process = subprocess.Popen(  cogdraw_command.split(),stderr= subprocess.PIPE,stdout=  subprocess.PIPE  )
	stdout,stderr = cogdraw_process.communicate()	
	
	
	
	
	

		