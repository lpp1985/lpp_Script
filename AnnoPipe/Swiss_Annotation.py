#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/12/6
"""
from Dependcy import *
from optparse import OptionParser
from Taxon_GI_Parse import *
import re

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
	end_list= glob.glob(out_put_path+'/*.xls')
	if end_list:
		sys.exit()

	if not os.path.exists( out_put_path ):
		os.makedirs( out_put_path )
	README = open(out_put_path+"Readme.txt",'w')
	README.write(
r"""
将所有的基因序列比对到Swwiss数据库,结果说明如下：

*.xls  详细的比对结果，用Excel打开。


"""	
	
	
	    )
	diamond_result = output_prefix+'_SwissAlignment.xls'
	error = RunDiamond(options.input,options.evalue, blast_type,"swissprot",diamond_result)
	if error:
		print( colored("%s 's Swiss process in Diamond of Nr is error!!","red") )
		print(colored( error,"blue"  ))
		print(  "##############################################"   )

		sys.exit()
	



