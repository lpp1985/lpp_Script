#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/11/14
"""

import os,sys
from os.path import abspath
from optparse import OptionParser
from lpp import *

usage = '''usage: python2.7 %prog -i input_path -t [The type you want]'''
parser = OptionParser(usage =usage ) 
parser.add_option("-o", "--Output", action="store",
                  dest="OutputPrefix",

                  help="OutputPrefix")

parser.add_option("-g", "--Path", action="store",
                  dest="Detail",
                  help="Pathway Detail File")

parser.add_option("-d", "--Diff", action="store",
                  dest="Diff",

                  help="Gene Different File")
parser.add_option("-a", "--Annotation", action="store",
                  dest="Anno",

                  help="Gene Annotation File")


parser.add_option("-e", "--Enrich", action="store",
                  dest="Enrich",

                  help="Gene Enrichment File")

if __name__ == '__main__':
    (options, args) = parser.parse_args()
    all_diff_gene = pd.read_table(options.Diff)
    annotation_data = pd.read_table( options.Anno  )
    all_diff_anno = annotation_data[ annotation_data["Name"].isin( all_diff_gene["id"]  )  ]
    
    ALL_PATHWAY = open( "%s"%(os.getpid()) ,'w')
    ALL_PATHWAY.write("Name\tPathwayID\tPathway\n")
    DETAIL = open(options.Detail,'rU')
    DETAIL.next()
    for line in DETAIL:
        line_l = line.split("\t")
        if not line_l[2]:
            continue
        name_list = line_l[2].split("||")
        for each_name in name_list:
            p_namelist = re.findall("(map\d{5})", each_name)
            for p_name in p_namelist:
                ALL_PATHWAY.write(line_l[0]+'\t'+p_name+'\t'+each_name+'\n')
            
    ALL_PATHWAY.close()
    all_pathway_data = pd.read_table( ALL_PATHWAY.name ).drop_duplicates()
    all_enriched_pathway  = pd.read_table( options.Enrich ).drop_duplicates()
    all_enriched_genePathway = all_pathway_data[  all_pathway_data["PathwayID"].isin (all_enriched_pathway["ID"])  ]
    enrich_go_annotation = pd.merge(
        left=all_enriched_genePathway, 
        right = all_diff_anno,
        left_on = "Name",
        right_on = "Name",
        how="inner"
    
    )    
    os.remove(ALL_PATHWAY.name)
    enrich_go_annotation = enrich_go_annotation.sort(["PathwayID"])
    enrich_go_annotation.to_csv(  options.OutputPrefix+".Annotation.tsv",sep='\t',index= False  )
    
    