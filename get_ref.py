from lpp import *
RAW = fasta_check(open(sys.argv[1],'rU'))
END = open(sys.argv[2],'w')
NEW = open(sys.argv[4],'w')
TAB = open(sys.argv[3],'w')
all_data = Ddict()
top_data = {}
all_new = {}
for t,s in RAW:
	gene = re.search("gene=(\S+)",t).group(1)
	if "ref_gene_id" not in t:
		all_new [gene] = ""

	rna_name = re.search("^\>(\S+)",t).group(1)
	all_data[gene][rna_name]=""
	if gene not in top_data:
		top_data[gene]=[len(s),rna_name,'>'+rna_name+'\n'+s]
	else:
		if len(s) >top_data[gene][0]:
			top_data[gene]=[len(s),rna_name,rna_name+'\n'+s]
				
				
for key in top_data:
	END.write('>'+key+'.'+top_data[key][-1])
	TAB.write( key+'\t'+ top_data[key][1]+'\t'+'; '.join( all_data[key] )   )
	if key in all_new:
		TAB.write("\tNew\n")
		NEW.write('>'+top_data[key][-1])
	else:
		TAB.write("\t-\n")
		
		