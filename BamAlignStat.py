#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/7/14
"""

from lpp import * 

all_align = Ddict()

if __name__ == '__main__':
	for a, b, c in os.walk("./"):
		for f in c:
			
			
			sample = f.split(".")[0]
			if f.endswith( "align.tsv" ):
				title = []
				for line in open(a + '/' + f):
					if "Stats" in line or ":" not in line:
						
						continue
					cate, value = re.split("\:\s+", line.strip())
					value = value.replace("\t", '')

					all_align[sample][cate] = value
					title.append(cate)
					
	END = open("TotalAlign.tsv", 'w')
	END.write("Sample\t" + "\t".join(title) + '\n' )
	for sample in all_align:
		END.write(sample)
		for cate in title:
			END.write("\t" + all_align[sample][cate])
		END.write("\n")
