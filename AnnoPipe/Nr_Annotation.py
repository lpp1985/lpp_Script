#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/12/6
"""
from Dependcy import *
from optparse import OptionParser
from Taxon_GI_Parse import *
import re

if __name__=="__main__":
	usage = '''usage: python2.7 %prog'''
	parser = OptionParser(usage =usage ) 
	parser.add_option("-i", "--INPUT", action="store", 
	                  dest="input", 
	                  help="input file")

	parser.add_option("-o", "--end", action="store", 
	                  dest="output_prefix", 
	                  help="output_prefix")

	parser.add_option("-e", "--evalue", action="store", 
	                  dest="evalue", 
	                  help="evalue cutoff")





	(options, args) = parser.parse_args() 
	FASTA = fasta_check(open(options.input,'rU'))
	sequence = FASTA.next()[-1]
	blast_type = Nul_or_Protein(sequence)
	output_prefix = os.path.abspath(  options.output_prefix )
	out_put_path = os.path.split(output_prefix)[0]+'/'
	end_list= glob.glob(out_put_path+'/*.xls')
	if end_list:
		sys.exit()

	if not os.path.exists( out_put_path ):
		os.makedirs( out_put_path )
	README = open(out_put_path+"Readme.txt",'w')
	README.write(
r"""
将所有的基因序列比对到Nr数据库,结果说明如下：

*.xls  详细的比对结果，用Excel打开。


"""	
	
	
	    )
	diamond_result = output_prefix+'_NrAlignment.tsv'
	error = RunDiamond(options.input,options.evalue, blast_type,"Nr",diamond_result)
	if error:
		print( colored("%s 's Nr process in Diamond of Nr is error!!","red") )
		print(colored( error,"blue"  ))
		print(  "##############################################"   )

		sys.exit()
	os.rename(output_prefix+'_NrAlignment.tsv',output_prefix+'_NrAlignment.xls')
	os.system("Nr_Taxon.py -i %s -o %s"%( output_prefix+'_NrAlignment.xls', output_prefix  )  )
	nr_data = pd.read_table(output_prefix+'_NrAlignment.xls')
	GENE_TAXON =  open( output_prefix+"_Taxon.txt",'w' )
	GENE_STATS =  open( output_prefix+"_TaxonStats.txt",'w' )
	taxon_stat_hash = Ddict()
	for i in xrange(0,len(nr_data)):
		gi = re.search("gi\|(\d+)",nr_data.loc[i,"Nr_Hit"])
		if gi:
			gi = gi.group(1)
			taxon_gi_sql = Taxon_GI.select(Taxon_GI.q.GI==int(gi) )   
			if taxon_gi_sql.count():
				taxon_gi_sql = taxon_gi_sql[0]
				taxon_id = taxon_gi_sql.Taxon
				taxon_name_sql = TaxonName.select(TaxonName.q.Taxon==taxon_id)   
				
				taxon_name = taxon_name_sql[0].Name
				
				GENE_TAXON.write( nr_data.loc[i,"Name"] +'\t'+taxon_name+'\n'  )

				taxon_stat_hash[taxon_name][ nr_data.loc[i,"Name"] ]=""
	for key in sorted( taxon_stat_hash,key= lambda x: len( taxon_stat_hash[x]  )   )[::-1]:
		
		GENE_STATS.write(   key+'\t%s'%(  len( taxon_stat_hash[key]  )  ) +'\n'  )



