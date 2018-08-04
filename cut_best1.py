#!/usr/bin/python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/4/28
from lpp import *
from optparse import OptionParser

if __name__=="__main__":
	usage = '''usage: python2.7 %prog'''
	parser = OptionParser(usage =usage ) 
	parser.add_option("-i", "--INPUT", action="store", 
                      dest="input", 
                      help="input file")

	parser.add_option("-o", "--end", action="store", 
                      dest="output", 
                      help="output")

	parser.add_option("-f", "--filter", action="store_true", 
                      dest="Filter",

	                  default=False,
                      help="For blast_Filter?")
	
	(options, args) = parser.parse_args() 
	RAW = open( options.input )
	END = open(options.output,'w')
	title = RAW.next()
	Filter = options.Filter
	if Filter:
		title_list = title.rstrip().split("\t")
		title_list =title_list[2:-3]
		title_list =[ "Nt_"+x.split('_',1)[-1] for x in title_list  ]
		title_list[0]="Name"
		END.write( "\t".join(title_list)+'\n'  )
	else:
		END.write(title)
	
	already = Ddict()
	#END.write(title)
	for line in RAW:
		line_l = line.rstrip().split('\t')

		line_l_new = line_l[2:-3]
		if line_l_new[-4]=='1':
			line_l_new[-4]='+'
		else:
			line_l_new[-4]='-'
		line_l_new[0] = line_l_new[0].split()[0]
		already[line_l[2]][float(line_l[10])] = '\t'.join(line_l_new)
		
	for key in already:
		key2 = sorted(already[key])[-1]
		END.write(already[key][key2]+'\n')
