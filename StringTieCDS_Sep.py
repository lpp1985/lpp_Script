from lpp import *
RAW = fasta_check(open(sys.argv[1],'rU'))
END = open(sys.argv[2],'w')
NEW = open(sys.argv[4],'w')
TAB = open(sys.argv[3],'w')
all_data = Ddict()
top_data = {}
all_new = {}
for t,s in RAW:
	gene = re.search("gene=(\S+)",t)
	ref_gene_id = re.search( "ref_gene_id=(\S+)",t   )
	if gene:
		gene = re.search("gene=(\S+)",t).group(1)
	else:
		gene = t[1:].split()[0]
	if ref_gene_id:
		gene = ref_gene_id.group(1)
	if "ref_gene_id" not in t:
		all_new [gene] = ""

	rna_name = re.search("^\>(\S+)",t).group(1)
	all_data[gene][rna_name]=""
	if gene not in top_data:
		top_data[gene]=[len(s),rna_name,s]
	else:
		if len(s) >top_data[gene][0]:
			top_data[gene]=[len(s),rna_name,s]
				
				
for key in top_data:
	END.write('>'+key+'\n'+top_data[key][-1])
	TAB.write( key+'\t'+ top_data[key][1]+'\t'+'; '.join( all_data[key] )   )
	if key in all_new:
		TAB.write("\tNew\n")
		NEW.write('>'+key+'\n'+top_data[key][-1])
	else:
		TAB.write("\t-\n")
		
		