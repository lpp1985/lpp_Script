#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/8/21
"""

from lpp import * 
from optparse import OptionParser


if __name__ == '__main__':
	usage='''usage: python %prog [options]

        Repeat db annotation by Repbase'''
	parser = OptionParser(usage =usage )
	parser.add_option(
	    "-i", "--Input", action="store",
	    dest="fasta",
	    type='string',
	    help="RAW sequence")
	
	parser.add_option(
	    "-d", "--Database", action="store",
	    dest="Database",
	    type='string',
	    default = "/pub/SOFTWARE/Other/RepeatMasker/Libraries/RepeatMasker.lib", 
	    help="Database")	

	parser.add_option(
	    "-o", "--Anno", action="store",
	    dest="Anno",
	    type='string',
	    help="RepeatMasker gff file")
	(options, args) = parser.parse_args()
	os.system("blastall -v 1 -b 1 -K 1  -p blastn -i %s  -d  %s  -a 64 -m 8 | cut -f 1,2|uniq  > db.list" % (options.fasta, options.Database) )
	anno_hash = {}
	for line in open("db.list"):
		line_l = line.split("\t")
		annotation = line.strip().split("#")[-1]
		anno_hash[line_l[0]] = annotation
		
	END = open(options.Anno, 'w')
	for t, s in  fasta_check(open( options.fasta, 'rU')):
		t =  t[1:-1].split()[0]
		if t in anno_hash:
			
			annotation =  anno_hash[t ]
		else:
			annotation = "Unknown"
		t = '>' + t + '#' + annotation + '\n'
		END.write(t + s)
		
		
		
		