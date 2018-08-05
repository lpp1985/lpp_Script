#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/12/30
"""
import sys,shlex,os,subprocess

from lpp import *
from optparse import OptionParser
from GO_obo_parse import *

usage = '''usage: python2.7 %prog -i input_path -t [The type you want]'''
parser = OptionParser(usage =usage ) 
parser.add_option("-i", "--INPUT", action="store", 
                  dest="input",

                  help="Input File")



parser.add_option("-o", "--end", action="store", 
                  dest="output", 
                  help="OUTPUT Data")
if __name__ == '__main__':
	(options, args) = parser.parse_args() 
	id_go = Ddict()
	for line in open(options.input,'rU'):
		line_l = line.strip().split("\t")
		for  e_f in line_l[1:]:
			id_go[line_l[0]] [e_f] = ""
	END = open(options.output,'w')
	END.write('Name\t'+"GO-BiologicalProcess"+'\t'+"GO-MolecularFunction"+"\tGo-CellularComponent\n")
	end_order = ["biological_process","molecular_function","cellular_component"]
	for key in id_go:
		all_go = id_go[key]
		mapp_result = Ddict()
		for each_go in all_go:
			component = GO_COMPONENT.select(GO_COMPONENT.q.Go==each_go)
			if not component.count():
				continue
			component = GO_COMPONENT.select(GO_COMPONENT.q.Go==each_go)[0].Compent
			
			mapp_result[ component ] [ each_go ]=""
	
		result = "\t".join([ '; '.join([ each_go+'//%s'%(  GO_DEF.select( GO_DEF.q.Go== each_go )[0].Def   )  for each_go in   mapp_result[x]  ] )     for x in end_order    ])
		END.write(key+"\t"+result+'\n')
	
