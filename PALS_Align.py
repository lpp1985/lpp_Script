#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/8/17
"""

from lpp import * 
from itertools import combinations
import os

if __name__ == '__main__':
	RAW = fasta_check( open(sys.argv[1], 'rU'))
	COMMAND = open( "command.list", 'w')
	all_hash = {}
	all_file = []
	i = 0
	for t, s in RAW:
		s1 = re.sub("\s+","",s)
		if len(s1)<100000:
			continue
		i += 1
		END = open( os .path.abspath("%s.cache" % (i)), 'w')
		END.write(t + s)
		END.close()
		all_file.append(END.name)
		COMMAND.write( "pals -self %s -out %s_self_cache.xml\n" % (END.name, i))
	i = 0
	for comb in combinations(all_file, 2):
		i += 1
		COMMAND.write("pals -query %s -target %s -out %s_target.xml\n" % (comb[0], comb[1], i))
		
	COMMAND.close()
	qstat = os.popen("which qstat")
	if not qstat:
		
		os.system("cat %s |parallel -j 60 " % (COMMAND.name ))
		all_file = glob.glob("*.xml")
		END = open(sys.argv[2], 'w')
		for e_f in all_file:
			
			END.write( open(e_f).read())
			os.remove(e_f)
	
		all_cache = glob.glob("*.cache")
		for e_f in all_cache:
			os.remove(e_f)
	else:
		os.system("Pals.nf --command %s --out %s" % (COMMAND.name,sys.argv[2] ))
