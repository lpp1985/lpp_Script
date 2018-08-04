#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/7/19
"""

from lpp import * 



if __name__ == '__main__':
	RAW = open(sys.argv[1], 'rU')
	SNP = open(sys.argv[2], 'w')
	FRAMESHIFT = open( sys.argv[3], 'w')
	title =  RAW.next()
	SNP.write(title)
	FRAMESHIFT.write( title)
	RAW = open(sys.argv[1], 'rU')
	data = RAW.read()
	data_l = re.split("\n(?=\S+)", data)
	for data_b in data_l:
		gene_name = data_b.split()[0]
		if "HIGH" in data_b:
			
			line_list = data_b.split("\n")
			for each_l in line_list:
				
				if "HIGH" in each_l:
					line_l = each_l.split("\t")
					line_l[0] = gene_name
					FRAMESHIFT.write("\t".join(line_l) + '\n')
		else:
			line_list = data_b.split("\n")
			for each_l in line_list:

				line_l = each_l.split("\t")
				line_l[0] = gene_name
				SNP.write("\t".join(line_l) + '\n')			
