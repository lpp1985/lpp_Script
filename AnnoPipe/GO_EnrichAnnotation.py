#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/3/16
"""
from lpp import *
from GO_obo_parse import *
usage = "python2.7 %prog [options]"
parser = OptionParser(usage =usage )
parser.add_option("-o", "--Output", action="store",
                  dest="OutputPrefix",

                  help="OutputPrefix")

parser.add_option("-g", "--GO", action="store",
                  dest="ALLGO",
                  help="AllGO Mapping File")

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
	
	#Extract All Differental Gene Annotation
	all_diff_gene = pd.read_table(options.Diff)
	all_Annotaion = pd.read_table( options.Anno )
	all_diff_anno = all_Annotaion[ all_Annotaion["Name"].isin( all_diff_gene["id"]  )  ]

	CACHE = open("%s.cache"%(os.getpid()) ,'w')
	RAW  = open(options.ALLGO)
	CACHE.write(RAW.next()[:-1]+'\tGODefination\n')
	
	for line in RAW:
		line_l =line.strip().split("\t")
		for key in line_l[1:]:
			define = GO_DEF.select( GO_DEF.q.Go== key )[0].Def
			CACHE.write(line_l[0]+'\t'+key+'\t'+define+'\n')
	CACHE.close()
	all_go = pd.read_table( CACHE.name )
	all_enrichgo = pd.read_table( options.Enrich )
	all_enrichgo_gene = all_go[ all_go["GOTerm"].isin( all_enrichgo["category"]  )  ]
	
	
	enrich_go_annotation = pd.merge(
	    left=all_enrichgo_gene, 
	    right = all_diff_anno,
	    left_on = "Gene",
	    right_on = "Name",
	    how="inner"
	    
	)
	del  enrich_go_annotation["Name"]
	enrich_go_annotation.to_csv(options.OutputPrefix+'.Annotation.tsv',sep = "\t",index = False)
	os.remove(CACHE.name)
	
