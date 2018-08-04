#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/9/26
"""

from lpp import * 
gene_dict = Ddict()


if __name__ == '__main__':
	END = open(sys.argv[2], 'w')
	all_data = open(sys.argv[1], 'rU').read()
	all_block = all_data.strip().split("\n\n")
	print(len(all_block))
	for e_b in all_block:
		data_l = e_b.split("\t")
		scaff = data_l[0]
		start, end = sorted((int(data_l[3]), int(data_l[4])))
		gene_dict[scaff][start] = e_b
		
		
	for scaff in sorted(gene_dict, key = lambda x: int(re.search("(\d+)", x).group(1))):
		for loc in sorted(gene_dict[scaff]):
			END.write( gene_dict[scaff][loc] + '\n\n\n')
		
