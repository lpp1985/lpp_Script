#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/1/11
"""
import sys
from lpp import *
from optparse import OptionParser


if __name__ == '__main__':
	usage='''usage: python %prog [options]

        Transform RepeatMasker gff result'''
	parser = OptionParser(usage =usage )
	parser.add_option(
	    "-i", "--Input", action="store",
	    dest="fasta",
	    type='string',
	    help="RAW sequence")

	parser.add_option(
	    "-g", "--gff", action="store",
	    dest="gff",
	    type='string',
	    help="RepeatMasker gff file")    

	parser.add_option(
	    "-s", "--seq", action="store",
	    dest="seq",
	    type='string',
	    help="Repeat Seuqence")

	(options, args) = parser.parse_args()

	GFF = open( options.gff,'rU'  )

	REPSEQ = open(options.seq, 'w')
	SEQ_ID = Ddict()
	i = 0
	for line in GFF:
		line_l = line.split("\t")
		if "Simple_repeat" in line_l[2] or "Low_complexity" in line_l[2] or "LTR_Repeat" in line_l[2]:
			continue
		i += 1
		seq_id = "Rep_%s" % (i )

		start = int(line_l[3])
		end = int(line_l[4])
		if start > end:
			end,start = start,end
		name = line_l[0].strip()
		SEQ_ID[name][seq_id] = [start, end]

	RAW = fasta_check(   open( options.fasta,'rU') )
	for t,s in RAW:
		name = t.split()[0][1:]
		s = re.sub("\s+","",s)
		for each_id, cont in SEQ_ID[name].items():
			REPSEQ.write(">%s\n%s\n" % (each_id, s[cont[0]:cont[1] ]) )
	REPSEQ.close()
	#os.system(" cdhit-est -c 0.8 -i  %s -o Cluster -M 0  -T 64" % (REPSEQ.name))
	
