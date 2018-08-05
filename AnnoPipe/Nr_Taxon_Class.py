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
parser.add_option("-n", "--Names", action="store", 
                  dest="Names", 
                  default = "",
                  help="Taxon Names")

(options, args) = parser.parse_args() 
NAME = open( options.Names, 'rU')
name_taxon = {}
for line in NAME:
    line_l = re.split("\s+\|\s+", line)
    name_taxon[ line_l[1] ] = int( line_l[0] )



nr_data = pd.read_table(options.Nr)

if not os.path.exists(path):
    os.makedirs(path)
GENE_TAXON =  open( options.OUT,'w' )

taxon_stat_hash = Ddict()
for i in xrange(0,len(nr_data)):
    taxon_name = re.search("\[([^\]]+)\]$",nr_data.loc[i,"Nr_Hit"])
    if taxon_name:
        taxon_name = taxon_name.group(1)

        taxon_id = name_taxon[ taxon_name]
        try:
            SuperKingdom = Taxon_Classification.select(Taxon_Classification.q.Taxon==taxon_id )[0].Class
            GENE_TAXON.write( nr_data.loc[i,"Name"] +'\t'+taxon_name+ '\t' +SuperKingdom + '\n'  )

            taxon_stat_hash[taxon_name][ nr_data.loc[i,"Name"] ]=""
        except:
            pass


    
