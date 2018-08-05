#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/3/7
"""
from lpp import *
RAW = open(sys.argv[1],'rU')
all_spieces = Ddict()
for line in RAW:
	line_l = line.split("\t")
	if line_l[2]!=line_l[0]:
		all_spieces[line_l[2]][line_l[0]]=""
def get_taxon_seed2( taxon_number ):
	all_end={}
	def creeper(taxon_number):
		all_end[ taxon_number ]=''
		if taxon_number in all_spieces:
			all_end[ taxon_number ]=''
			for key1 in all_spieces[ taxon_number ]:
				creeper(key1)
	creeper(taxon_number)
	return all_end
all_bac = get_taxon_seed2("2")
END = open("Taxon_class",'w')
for key in all_bac:
	END.write(key+'\tBacter\n')
all_Animal = get_taxon_seed2("33208")
all_Plants = get_taxon_seed2("3193")
for key in all_Animal:
	END.write(key+'\tAnimal\n')

for key in all_Plants:
	END.write(key+'\tPlants\n')
