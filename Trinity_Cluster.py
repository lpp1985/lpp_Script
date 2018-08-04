#!/usr/bin/python
from lpp import *
RAW = fasta_check(open(sys.argv[1],'rU'))
comp_data = Ddict()
for t,s in RAW:
	t = t.split()[0][1:]
	comp = re.search("(^.*c\d+_g\d+)",t).group(1)
	comp_data[comp][t]=s
SEQ = open( "Trinity_cluster.fa",'w'  )
LIST = open("Trinity_cluster.tsv",'w')
LENGTH = open("Trinity_clusterLength.tsv",'w')
for each_comp in comp_data:
	clusterid = sorted(  comp_data[each_comp],key=lambda x: len(comp_data[each_comp][x] )   )[-1]
	seq = comp_data[each_comp][clusterid]
	SEQ.write( '>'+clusterid+'\n'+ seq)
	LIST.write( each_comp+'\t'+'; '.join( comp_data[each_comp]  )+'\t'+clusterid+'\n'    )
	LENGTH.write( each_comp+'\t'+str( len(re.sub("\s+","",seq)) ) +'\n'    )
	
