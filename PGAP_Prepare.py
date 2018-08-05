#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2018/1/5
"""

from lpp import * 



if __name__ == '__main__':
	all_path = sys.argv[1].split("/")[0]
	
	ALLFAA = fasta_check( open(all_path + '/' + all_path + ".faa"))
	all_seq = {}
	FAA = open(all_path + '.pep', 'w' )
	for t, s in ALLFAA:
		name = t.split()[0][1:]
		all_seq[name] = ""
		FAA.write(t.split()[0] + '\n' + s)
		
	ALLNUC = fasta_check( open(all_path + '/' + all_path + ".ffn"))
	NUC = open(all_path + '.nuc', 'w' )
	for t, s in ALLNUC:
		name = t.split()[0][1:]
		if name in all_seq:
			
			NUC.write(t.split()[0] + '\n' + s)
	FUNCTION = open(all_path + ".function", 'w')
	for line in open(all_path + '/Annotation/Table/GeneFeature+Annotation.xlsx'):
		line_l = line.split("\t")
		if line_l[0] in all_seq:
			if "COG" in line_l[59]:
				cog = line_l[59] + line_l[61]
			else:
				cog = "-"
			FUNCTION.write(line_l[0] + '\t' + cog + '\t' + line_l[2] + '\n')