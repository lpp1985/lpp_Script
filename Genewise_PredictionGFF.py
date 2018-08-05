#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/10/20
"""

from lpp import *
from optparse import OptionParser



if __name__ == '__main__':
	usage = '''usage: python2.7 %prog [options] '''
	parser = OptionParser(usage =usage )    

	parser.add_option("-i", "--genewise", action="store",
                      dest="genewise",
                      help="Genese Result!!")

	parser.add_option("-o", "--output", action="store",
                      dest="output",
                      help="GFF Parse output!!")

	parser.add_option("-t", "--train", action="store",
                      dest="training",
                      help="GFF for Braker_new Training!!")
	parser.add_option("-p", "--protein", action="store",
	                  dest="protein",
	                  help="Protein sequence !!")	

	(options, args) = parser.parse_args()
	protein_length = {}
	for t, s in fasta_check( open( options.protein, 'rU')):
		s = re.sub("\s+", "", s)
		name = re.search("(STRG\.[^\:]+)", t).group(1)
		protein_length[name] = len(s)

	genewise = options.genewise
	output = options.output  
	END = open(output,'w')
	CACHE = open("Cache.gff3", 'w')
	TRAINING = open( options.training, 'w')
	RAW = re.split("//\nBits",open(genewise,'rU').read())
	j=0
	for e_b in RAW:
		if "STRG" not in e_b:

			continue
		mrna_name = re.search("(STRG\S+)", e_b).group(1)
		gene_name =  mrna_name.rsplit(".", 1)[0]
		need_data=  re.search("\n\S+\s+STRG\S+\s+(\d+)\s+(\d+)", e_b  )
		all_start = need_data.group(1)
		all_end =  need_data.group(2)
		if all_start !="1":
			continue
		if int(all_end) != protein_length[mrna_name]:
			continue
		data_b = e_b.split("\n//\n")
		raw_b = data_b[-1]

		all_b = re.split("\n(?=.+match)", raw_b)
		raw_score = 0 
		for e_all_b in all_b:

			if not e_all_b.strip():
				continue
			score = re.search("match\s+\d+\s+\d+\s+(\S+)",e_all_b ).group(1)
			if float(score) > raw_score:
				gff_b =  e_all_b
				raw_score = score
		cds_data = []
		exon_data = []


		for line in gff_b.split("\n"):
			if not line:
				continue

			gff_list = line.split("\t")
			if gff_list[2] == "intron":
				continue


			seq_name = gff_list[0]

			subjname,loc_append = seq_name.rsplit("__",1)
			loc_append  = int(loc_append)
			gff_list[0] = subjname
			gff_list[3] = int(gff_list[3]) + loc_append
			gff_list[4] = int(gff_list[4]) + loc_append


			gff_list[3],gff_list[4] = sorted([ gff_list[3],gff_list[4] ])
			gff_list[3] = str(gff_list[3])
			gff_list[4] =  str(gff_list[4])


			if gff_list[2] == "match":
				gff_list[2] = "mRNA"
				if gff_list[6]=="+":
					gff_list[4] = str( int(gff_list[4]) + 3)
				else:
					gff_list[3] = str( int(gff_list[3]) - 3)
				gff_list[-1] = "ID=%s;Parent=%s" % (mrna_name,gene_name)
				CACHE.write("\t".join(gff_list)+'\n')
			if gff_list[2] == "cds":
				gff_list[2] = "exon"
				exon_data.append("\t".join(gff_list[:-1]))
				gff_list[2] = "CDS"
				cds_data.append( "\t".join(gff_list[:-1]))

		cds_data = sorted( cds_data, key = lambda x: int( x.split("\t")[3] ) ) 
		exon_data =  sorted( exon_data, key = lambda x: int( x.split("\t")[3] ) )
		if gff_list[6] == "+":
			new_list =   cds_data[-1].split("\t")
			new_list[4] = str( int(new_list[4]) + 3)
			cds_data[-1] = "\t".join( new_list)
			new_list =   exon_data[-1].split("\t")
			new_list[4] = str( int(new_list[4]) + 3)
			exon_data[-1] = "\t".join( new_list)

		else:
			new_list =   cds_data[0].split("\t")
			new_list[3] = str( int(new_list[3]) - 3)
			if int(new_list[4])<0:
				continue
			cds_data[0] = "\t".join( new_list)
			new_list =   exon_data[0].split("\t")
			new_list[3] = str( int(new_list[3]) - 3)
			exon_data[0] = "\t".join( new_list)            
		cds_data = iter(cds_data )
		exon_data = iter(exon_data )
		i = 1
		for key in exon_data:
			CACHE.write(key + '\t' + "ID=%s.Exon%s;Parent=%s\n" % (mrna_name, i, mrna_name))

			CACHE.write(cds_data.next() + '\t' + "ID=%s.CDS%s;Parent=%s\n" % (mrna_name, i, mrna_name))
			i += 1
		CACHE.write("\n")
	CACHE.close()

	CACHE_RAW = open(CACHE.name, 'rU')
	cache_data = CACHE_RAW.read()
	cache_data_list = cache_data.split("\n\n")
	all_gene_dict = Ddict()
	CACHE2 = open( "Cache2.gff", 'w')

	for cache_block in  cache_data_list:
		cache_block =  cache_block.strip()
		if not cache_block:
			continue

		gene_name = re.search("(STRG\.\d+)", cache_block ).group(1)
		if gene_name not in all_gene_dict:
			TRAINING.write(cache_block + '\n\n')
		all_gene_dict[gene_name][cache_block] = ""


	for each_gene, gene_gff in all_gene_dict.items():
		all_loc = []
		for each_con in gene_gff:
			[ (start, end) ]= re.findall("mRNA\t(\d+)\t(\d+)", each_con)
			all_loc.append( int(start) )
			all_loc.append( int(end) )

		data = re.search("(.+\tmRNA\t.+)", each_con).group(1)    
		start = str(sorted(all_loc)[0])
		end =  str( sorted(all_loc)[-1] )
		data_list = data.split("\t")
		data_list[2] = "Gene"
		data_list[3] = start
		data_list[4] =  end
		data_list[5] = "0"
		data_list[-1] = "ID=%s;Name=%s" %(each_gene,each_gene )


		CACHE2.write("\t".join(data_list) + '\n')
		CACHE2.write( "\n".join(gene_gff) + '\n\n' )

	CACHE2.close()

	all_data = open(CACHE2.name, 'rU').read()
	all_data_list = all_data.split("\n\n")
	all_data_hash = Ddict()
	for e_b in all_data_list:
		if not e_b.strip():
			continue
		name_l = e_b.split("\t")
		scf = name_l[0]
		start = int(name_l[3])
		all_data_hash[scf][start][e_b] = ""


	for scf in sorted(all_data_hash):
		for loc in sorted(all_data_hash[scf]):
			for e_b in  all_data_hash[scf][loc]:
				END.write(e_b + '\n')
	os.remove(CACHE2.name)
	os.remove(CACHE.name)


