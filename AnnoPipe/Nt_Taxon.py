#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2016/3/9
"""
from Dependcy import *
from optparse import OptionParser
from Taxon_GI_Parse import *
import re
usage = '''usage: python2.7 %prog'''
parser = OptionParser(usage =usage ) 
parser.add_option("-i", "--Input", action="store", 
                  dest="Nr", 
                  default = "",
                  help="Nr XLS file")
parser.add_option("-o", "--Out", action="store", 
                  dest="OUT", 
                  default = "",
                  help="Output prefix")
(options, args) = parser.parse_args() 


nr_data = pd.read_table(options.Nr)
output_prefix = os.path.abspath( options.OUT  )
path = os.path.dirname(output_prefix)
if not os.path.exists(path):
    os.makedirs(path)
GENE_TAXON =  open( "%s_Taxon.txt"%(output_prefix),'w' )
GENE_STATS =  open( "%s__TaxonStats.txt"%(output_prefix),'w' )
GENE_STATS.write(  "Taxon\tNumber\tPercentage\n"  )
taxon_stat_hash = Ddict()
for i in xrange(0,len(nr_data)):
    gi = re.search("gi\|(\d+)",nr_data.loc[i,"Nt_Hit"])
    if gi:
        gi = gi.group(1)
        taxon_gi_sql = Taxon_GI.select(Taxon_GI.q.GI==int(gi) )   
        if taxon_gi_sql.count():
            taxon_gi_sql = taxon_gi_sql[0]
            taxon_id = taxon_gi_sql.Taxon
            taxon_name_sql = TaxonName.select(TaxonName.q.Taxon==taxon_id)   
            try:
                taxon_name = taxon_name_sql[0].Name
		SuperKingdom = Taxon_Classification.select(Taxon_Classification.q.Taxon==taxon_id )[0].Class
    
                GENE_TAXON.write( nr_data.loc[i,"Name"] +'\t'+taxon_name+'\t'+SuperKingdom+'\n'  )

                taxon_stat_hash[taxon_name][ nr_data.loc[i,"Name"] ]=""
            except:
                pass
total = 0
for taxon in taxon_stat_hash:
    total+=len(taxon_stat_hash[taxon])
    
for key in sorted( taxon_stat_hash,key= lambda x: len( taxon_stat_hash[x]  )   )[::-1]:

    GENE_STATS.write(   key+'\t%s'%(  len( taxon_stat_hash[key]  )  ) +'\t'+str(  1.0*len( taxon_stat_hash[key]  )/total   )+'\n'  )

    
