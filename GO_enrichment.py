#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/11/4
"""

from enrichment import *
from optparse import OptionParser
usage = '''usage: python %prog [options] 

It can automaticly do GO Mapping & Enrichment Analysis!!'''


parser = OptionParser(usage =usage ) 

parser.add_option("-d", "--DIFF", action="store", 
                  dest="diff",
                  type='string',
                  help="the Difference GO Mapping Result")

parser.add_option("-a", "--ALL", action="store", 
                  dest="all",
                  type='string',
                  help="the ALL GO Mapping Result")


parser.add_option("-o", "--OUT", action="store", 
                  dest="out",
                  type='string',
                  help="Result")
(options, args) = parser.parse_args() 
difference_mapping = options.diff

END = open(options.out,'w')
all_mapping = options.all
ALL = open(all_mapping,'rU')
ALL.next()
go_def = {}
all_anno = 0
all_diff = 0
go_all_number = {}
go_diff_number = {}
for line in ALL:
	if not line.startswith("GO"):
        	continue
	line_l = line.strip().split("\t")
	go_def[line_l[0]] =line_l[1] 
	go_all_number[line_l[0]] = int(line_l[2])
	all_anno+=int(line_l[2])
	go_diff_number[line_l[0]] = 0
DIFF = open(difference_mapping,'rU')
DIFF.next()
END.write("Go\tDef\tTotal\tDiff\tp_val\tq_val\n")

for line in DIFF:
	if not line.startswith("GO"):
		continue
	line_l = line.strip().split("\t")
	go_diff_number[line_l[0]] += int(line_l[2])
	all_diff+=int(line_l[2])
p_val_list = []
for key in sorted(go_diff_number):

	p_val_list.append(enrichment_analysis(all_diff, all_anno, go_all_number[key] , go_diff_number[key]))


q_val_iter = iter(fdr(p_val_list))
p_val_iter = iter(p_val_list)

for key in sorted(go_diff_number):
	q_val = q_val_iter.next()
	p_val = p_val_iter.next()

	if q_val<0.05:
		END.write(key+'\t'+go_def[key]+"\t%s\t%s"%(go_all_number[key],go_diff_number[key])+"\t%s\t%s\n"%(p_val,q_val))
		
