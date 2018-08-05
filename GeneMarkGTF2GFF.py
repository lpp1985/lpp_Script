#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/1/16
"""

import sys,re
from optparse import OptionParser


if __name__ == '__main__':
	usage='''usage: python %prog [options]

	    conver gtf to gff3 from genmarks result 1 '''
	parser = OptionParser(usage =usage )

	parser.add_option("-i", "--GTF", action="store",
	                  dest="gtf",
	                  type='string',
	                  help="gtf file")




	parser.add_option("-o", "--GFF", action="store",
	                  dest="gff",

	                  help="Result gff")
	(options, args) = parser.parse_args()
	RAW = open(options.gtf,'rU')
	END = open(options.gff,'w')
	gene_id = ""
	gene_hash = {}
	direction  = "+"
	chro = ""
	for line in RAW:
		line_l = line.split("\t")
		gene_id2 = re.search("gene_id\s+\"(\S+)\"",line).group(1)
		Product2 =  re.search("Product\s+\"([^\"]+)\"" ,line).group(1)
		transcript_id2 = re.search("transcript_id\s+\"(\S+)\"",line).group(1)
		direction2 = line_l[6]
		if gene_id=="":
			chro = line_l[0]
			gene_id = gene_id2
			direction = direction2
			Product = Product2
			transcript_id = transcript_id2
		if gene_id2 ==gene_id:
			if line_l[2]=="CDS":
				gene_hash[int(line_l[3])] = "\t".join(  line_l[:-1]  )
				gene_hash["end"] = line_l[4]
		else:
			if line_l[2]=="CDS":
				cds_number = len(gene_hash)-1
				gene_max = gene_hash["end"]
				if direction =='-':
					CDS_number = len(gene_hash)
				else:
						CDS_number=0				
				del gene_hash["end"]
				
				END.write( chro+'\tGeneMark.hmm\tgene\t' + str( sorted( gene_hash)[0] )+ '\t'+str( gene_max  ) +'\t'+direction+'.'   )
				END.write("\tID=%s; Product=\"%s\";\n"%(gene_id,Product))
				for each_coor in sorted(gene_hash):
					END.write(   gene_hash[ each_coor  ]  +'\t' )
					if direction=='+':
						CDS_number+=1
					else:
						CDS_number-=1
					END.write(   "\tID=%sCDS%s; Transcript=%s; Parent=%s; Product=\"%s\";\n"%(gene_id,CDS_number,transcript_id,gene_id,Product)   )
				
				gene_hash = {}
				gene_hash[int(line_l[3])] = "\t".join(  line_l[:-1]  )
				gene_hash["end"] = line_l[4]
				gene_id = gene_id2
				direction = direction2
				Product = Product2	
				chro = line_l[0]
				cds_number = len(gene_hash)-1
	gene_max = gene_hash["end"]
	if direction =='-':
		CDS_number = len(gene_hash)
	else:
		CDS_number=0				
	del gene_hash["end"]

	END.write( chro+'\tGeneMark.hmm\tgene\t' + str( sorted( gene_hash)[0] )+ '\t'+str( gene_max  ) +'\t'+direction+'.'   )
	END.write("\tID=%s; Product=\"%s\";\n"%(gene_id,Product))
	for each_coor in sorted(gene_hash):
		END.write(   gene_hash[ each_coor  ]  +'\t' )
		if direction=='+':
			CDS_number+=1
		else:
			CDS_number-=1
		END.write(   "\tID=%sCDS%s; Transcript=%s; Parent=%s; Product=\"%s\";\n"%(gene_id,CDS_number,transcript_id,gene_id,Product)   )
	