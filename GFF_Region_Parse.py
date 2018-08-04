#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/7/17
"""

from lpp import * 



if __name__ == '__main__':
	RAW = open(sys.argv[1], 'rU')
	data = RAW.read()
	total_protein_loc = Ddict()
	o_proteinid = ""
	o_scaff = ""
	RAW = open(sys.argv[1], 'rU')
	if "locus_tag" in data and "YAL068C" in data:
		for line in RAW:
			if "\tgene\t" in line:
				line_l = line.split("\t")
				protein_id = re.search("locus_tag\=([^\;\s]+)", line).group(1)
				[start, end] = sorted( [int(line_l[3]), int(line_l[4])] )
				total_protein_loc[ line_l[0] ][start][end] = protein_id
				
		
	elif "protein_id"  in data:
		location_list = []
		for line in RAW:
			line_l = line.split("\t")
			if "\tCDS\t" in line:
				scaff = line_l[0]
				protein_id= re.search("protein_id\=([^\;\s]+)", line)
				if not protein_id:
                                      
					continue
				protein_id = protein_id.group(1)
				if o_scaff == "":
					
					o_scaff = scaff
                 
				if o_proteinid != protein_id and o_proteinid != "":
					location_list = sorted(location_list)
					total_protein_loc[o_scaff][location_list[0]][location_list[-1]] = o_proteinid
			
					o_proteinid = protein_id
					line_l = line.split("\t")
					o_scaff = scaff
					location_list = [ int(line_l[3]), int( line_l[4] )]
				else:
					location_list.extend ([ int(line_l[3]), int( line_l[4] )])
					if o_proteinid == "":
						o_proteinid = protein_id
		location_list = sorted(location_list)
		total_protein_loc[o_scaff][location_list[0]][location_list[-1]] = o_proteinid
	else:
		for line in RAW:
			if "\tgene\t" in line:
				line_l = line.split("\t")
				gene_id = re.search("ID\=gene_(\d+)", line).group(1)
				protein_id =  "CDS_" + gene_id
				[start, end] = sorted( [int(line_l[3]), int(line_l[4])] )
				total_protein_loc[ line_l[0] ][start][end] = protein_id
END = open( "All_CDS.txt", 'w')
for chrom in total_protein_loc:
	for start in sorted( total_protein_loc[ chrom]):
		
		for end in  sorted( total_protein_loc[ chrom][start]):
			END.write( chrom + '\t%s\t%s\t%s\n' % (start, end,total_protein_loc[ chrom][start][end] ))
