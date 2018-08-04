#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/8/22
"""

from lpp import * 
from optparse import OptionParser


if __name__ == '__main__':
	usage='''usage: python %prog [options]

        Transform RepeatMasker gff result'''
	parser = OptionParser(usage =usage )
	parser.add_option(
	    "-i", "--Input", action="store",
	    dest="input",
	    type='string',
	    help="repeat masker out file")

	parser.add_option(
	    "-o", "--output", action="store",
	    dest="output",

	    help="output ")
	parser.add_option(
	    "-s", "--all_seq", action="store",
	    dest="seq",

	    help="all sequence in fasta ")	
	(options, args) = parser.parse_args()
	RAW = open( options.input)
	all_result = Ddict()
	one_cate = {}
	total_length = 0
	for line in RAW:
		if not re.match("^\s*\d+", line):
			
			continue
		line_l = line.split()
	
		length = abs(int(line_l[6]) - int(line_l[5])) + 1
		total_length += length
		category = line_l[10].replace("?", "").split("/")
		if len(category) == 2:
			if all_result[ category[0]][category[1]]:
				all_result[ category[0]][category[1]] += length
			else:
				all_result[ category[0]][category[1]] = length
		else:
			if category[0] not in  one_cate:
				one_cate[category[0]  ] = length
			else:
				one_cate[category[0]  ] += length
	all_number = 0
	for t, s in fasta_check( open(options.seq, 'rU')):
		all_number += len(re.sub("\s+", "", s))
	END = open(options.output, 'w')
	END.write("Category\tSubCategory\tLength\tPerc(%)\n")
	all_del = []
	for each_category in sorted(one_cate):
		if each_category in all_result:
			all_result[each_category]["Unknown"] =  one_cate[each_category ]
			all_del.append( each_category)
	
	
	
	for each_category in all_result:
		total = 0
		for each_sub in all_result[each_category]:
			total +=  all_result[each_category][ each_sub]
			perc = total * 100.0 / all_number
		END.write(each_category + '\ttotal\t%s\t%.2f\n' % (total, perc) )
		
		for each_sub in sorted(all_result[each_category]):
			perc =  all_result[each_category][ each_sub] * 100.0 / all_number
			END.write(  '\t%s\t%s\t%.2f\n' % (each_sub, all_result[each_category][ each_sub], perc) )
			
	for each_category in sorted(one_cate):
		if each_category in  all_del:
			continue
		number =  one_cate[each_category ]
		perc = number * 100.0 / all_number
		END.write(each_category + '\t\t%s\t%.2f\n' % (number, perc) )
	print(total_length)