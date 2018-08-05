#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/9/6
"""

from lpp import * 



if __name__ == '__main__':
	all_f =  glob.glob("*.stats")
	END = open("Assembly.xls", 'w')
	END.write( "Sample\tScaff.No\tN50\tLongest\tN90\tSum.Base\n")
	for e_f in all_f:
		name = e_f.split(".")[0]
		data = open(e_f).read()
		all_data = re.search("\n(\d+)\s+scaffolds from \d+ contigs sum up (\d+)bp, with average length", data)
		scff_no = all_data.group(1)
		sumbase = all_data.group(2)
		all_data = re.search("the longest is (\d+)bp,scaffold N50 is (\d+) bp, scaffold N90 is (\d+)\s+bp", data)
		longest = all_data.group(1)
		n50 = all_data.group(2)
		n90 = all_data.group(3)
		END.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (name, scff_no, n50, longest, n90, sumbase) )
