#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/1/13
"""

from lpp import *




if __name__ == '__main__':
	all_annotation = File_dict( open(sys.argv[1],'rU') ).read(1,2)
	GTF = open(  sys.argv[2]  ,'rU')
	END = open(  sys.argv[3],'w' )
	for line in GTF:
		name =  re.findall("gene_id \"(\S+)\"",line)
		if name:
			name = name[0]
			if name in all_annotation:
				function = all_annotation[ name ]
			else:
				function = "Hypothetical Protein"
			line = line.strip()+'Product \"%s\";\n'%(function)
			
		END.write(line)
