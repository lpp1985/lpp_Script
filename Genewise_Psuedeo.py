#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/1/13
"""

from lpp import *
from optparse import OptionParser
import os

if __name__ == '__main__':
	usage='''usage: python %prog [options]

	    multiproecssing blast '''
	parser = OptionParser(usage =usage )

	parser.add_option("-p", "--Protein", action="store",
		              dest="protein",
		              type='string',
		              help="protein file")

	parser.add_option("-n", "--Nucl", action="store",
		              dest="nucleotide",
		              type='string',
		              help="nucleotide file")


	parser.add_option("-o", "--Oputpu", action="store",
		              dest="out",

		              help="Result gff")		

	parser.add_option(
	    "-d","--rev",
	    action="store_true",
	    default=False,
	    dest="rev"
	                  )



	(options, args) = parser.parse_args()
	protein = options.protein
	PRO = open(protein)
	GFF = open(options.out,'w')
	rev = options.rev
	line  = PRO.next()
	anno =line [1:].strip().split(" ",1)[-1].strip()
	
	seq_name = line[1:].split()[0]
	nucl = options.nucleotide
	pu_number =  os.path.basename(nucl) .split(".nuc")[0].split("_")[-1]
	commandline = " genewise %s  %s -silent -quiet  "%( protein,nucl  )  
	if rev :
		commandline +="-trev"
	data = os.popen(commandline ).read()
	data = data.split("See WWW help for more info")[-1]
	gene_id  ="P__"+ seq_name+"__%s"%(pu_number)
	if "!"   in data or "X"  in data:
		commandline = " genewise %s  %s -silent -quiet -gff  "%( protein,nucl  )  
		if rev:
			commandline+=" -trev "
		GFF_RAW = os.popen(commandline  )
		exon_num=0
		intron_num = 0
		for line in GFF_RAW:
			if "//" in line:
				continue
			
			line_l = line.strip().split("\t")
			location = int(line_l[0].split("__")[-1])
			line_l[0] = line_l[0].rsplit("__")[0]
			line_l[3] = str(int(line_l[3])+location)
			line_l[4] = str(int(line_l[4])+location)
			if line_l[2] =="match":
				line_l[2] ="Pseudogene"
				line_l[-1]= "ID=%s;Name=%s;Product=\"%s[Pseudogene]\""%(gene_id,gene_id,anno)
			elif line_l[2]=="cds":
				exon_num+=1
				line_l[2] ="pseudogenic_exon"
				line_l[-1]= "ID=%sExon%s;Parent=%s;Product=\"%s[Pseudogene]\""%(gene_id,exon_num,gene_id,anno)
				
			else:
				intron_num+1
				line_l[2] ="pseudogenic_intron"
				line_l[-1]= "ID=%sIntron%s;Parent=%s"%(gene_id,intron_num,gene_id)
			GFF.write("\t".join(line_l)+'\n')
