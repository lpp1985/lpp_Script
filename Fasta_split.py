#!/usr/bin/env python
#coding:utf-8
"""
  Author:  LPP --<Lpp1985@hotmail.com>
  Purpose: 
  Created: 2015/1/3
"""
import os,sys

sys.path.append(os.path.split(__file__)[0]+'/../Lib/')
from lpp import *
from optparse import OptionParser

usage = "python2.7 %prog [options]"
parser = OptionParser(usage =usage )
parser.add_option("-i", "--Input", action="store",
                  dest="Input",

                  help="Input Fasta")
parser.add_option("-o", "--Output", action="store",
                  dest="output",

                  help="OutputPath")

parser.add_option("-t", "--Threshold", action="store",
                  dest="Threshold",
                  type="int",
                  help="Length Threshold!!!")

parser.add_option("-n", "--Name", action="store",
                  dest="Strain",
                  default="",
                  help="Strain Name!!!")

if __name__ == '__main__':
	(options, args) = parser.parse_args()
	Input = options.Input
	OUTPUT = options.output	
	Threshold = options.Threshold
	Strain = options.Strain
	g=0
	p=0
	all_seq = {}
	for t,s in fasta_check(open( options.Input,'rU') ):
		all_seq[s]=[t,s]
	for se in sorted( all_seq,key = lambda x: len(x)   )[::-1]:
		t,s  = all_seq[se]
		s2 = re.sub("\s+","",s)
		length = len(s2)
		status = re.search("(\w+)\n",t).group(1)
		if length>Threshold:
			g+=1
			Name=Strain+"_Genome%s"%(g)
		else:
			p+=1
			Name=Strain+"_Plasmid%s"%(p)
		OUT = open(OUTPUT+'/%s.fasta'%(Name),'w')
		OUT.write('>%s Length=%s %s\n%s'%(
	        Name,
	        length,
	        status,
	        s
	    )
	            )
		
			
