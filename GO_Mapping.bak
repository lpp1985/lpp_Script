#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/27
from lpp import *
from optparse import OptionParser
import operator as op
from enrichment import *
all_geneid_id = Ddict()

diff_geneid_id = Ddict()



go_id = Ddict()

go_def = {}

total_number_check = Ddict()
all_mapped_gene = {}
diff_mapped_gene = {}


'This function is call the son of each MF GO'
def recall( leaf, number_check ):
	global root_static,need_check,go_id
	for each_son in  son_father[ leaf ]:
		if each_son not in all_stop:

			recall( each_son, number_check )
		else:
			root_static[ each_son ][ need_check ] = ''
			for each_id in go_id[need_check]:
				number_check[ each_son ][ each_id ] = ''

'''Compute the combination value'''
usage = '''usage: python %prog [options] 

It can automaticly do GO Mapping & Enrichment Analysis!!'''


parser = OptionParser(usage =usage ) 

parser.add_option("-d", "--DIFF", action="store", 
                  dest="diff",
                  type='string',
                  help="the Difference GO Mapping Result")

parser.add_option("-a", "--ALL", action="store", 
                  dest="all",
                  type='string',
                  help="the ALL GO Mapping Result")
parser.add_option("-o", "--OUT", action="store", 
                  dest="output",
                  type='string',
                  help="the GO Mapping Result")
(options, args) = parser.parse_args() 
difference_mapping = options.diff
output_result = options.output

all_mapping = options.all



DIFF = open( difference_mapping,'rU'  )

DIFF.next() # delete the Header


ALL = open(  all_mapping,'rU'   )

ALL.next(  ) # delete the Header
all_diff = {}
all_go = {}
for line in DIFF:

	line_l = line.split('\t')

	name = line_l[2]

	gi = re.search( '(?:sp|tr)\|(\w+)' ,line_l[5] ).group(1)

	output = '-\t-\t->\t'+name+'\t'+gi


	diff_geneid_id[ gi ][ name ] = output
	all_diff[name] = ''
for line in ALL:

	line_l = line.split('\t')

	name = line_l[2]

	gi = re.search( '(?:sp|tr)\|(\w+)' ,line_l[5] ).group(1)
	all_go[name] = ''
	output = '-\t-\t->\t'+name+'\t'+gi


	all_geneid_id[ gi ][ name ] = output



ALL_GO = open( '/home/lpp/Database/GO/gene_association.goa_uniprot','rU'  )
ALL_GO.next()

DEF = open( '/home/lpp/Database/GO/NAME_DEF.list','rU'  )
change_go = File_Ddict( open(  '/home/lpp/Database/GO/relationship.alter','rU'  )  ).read(1,2)


for line in DEF:
	line_l = line.split('\t')
	go_def[ line_l[0] ] = line_l[1]

total_geneid_id = all_geneid_id
# combine the two hash
id_go = Ddict()
for line in ALL_GO:
	line_l = line.split('\t')
	if line_l[0] !='UniProtKB':
		continue
	upr_id = line_l[1]
	if upr_id in total_geneid_id:
		if upr_id in diff_geneid_id:

			for each_name in diff_geneid_id[  upr_id  ]:
				diff_mapped_gene[each_name] = ""
		for each_name in total_geneid_id[  upr_id  ]:
			all_mapped_gene[each_name] = ""

			if line_l[4] not in change_go:
				go_id[ line_l[4] ][ each_name ] = total_geneid_id[  upr_id  ][ each_name  ]+'\t'+line_l[4]+'\t'+go_def[  line_l[4] ]
				id_go[ each_name ][ line_l[4] ] = ''

			else:
				for each_altered in change_go[ line_l[4] ]:
					go_id[ each_altered ][ each_name ] = diff_geneid_id[  name  ][ each_name  ]+'\t'+each_altered+'\t'+go_def[  each_altered ]
					id_go[ each_name ][ each_altered ]


son_father = File_Ddict( open('/home/lpp/Database/GO/relationship.son','rU' ) ).read(2,1)
all_stop = File_Ddict(  open('/home/lpp/Database/GO/ROOT.root' ,'rU') ).read(1)
root_static = Ddict()



for each_go in go_id:
	need_check = each_go
	each_son = recall( each_go, total_geneid_id )

END = open( '%s.GO-mapping.list'%(output_result),'w'  )
END.write( 'ROOT\tROOT_DEF\tNumber\tID\tUnirpotKB\tGO\tGO_DEF\n'  )

number_all = len(all_mapped_gene)
number_diff = len(diff_mapped_gene)

p_value_list = []
enrich_cache = []


for each_go in reversed( sorted(root_static,key= lambda x: len( root_static[ x  ]) ) ):
	all_spec_go = {}
	diff_sepc_go = {}	
	for each_son in root_static[ each_go ]:
		for each_id in go_id[ each_son ]:
			all_spec_go[each_id] = ''	
			if each_id in all_diff:
				diff_sepc_go[each_id] = ""		
	END.write( each_go+'\t'+ go_def[  each_go ] +'\t%s\n'%( len(all_spec_go) )  ) 
	

	
	
	
	for each_son in root_static[ each_go ]:
		for each_id in go_id[ each_son ]:
			END.write(  go_id[ each_son ][  each_id ]+'\n' )
				
	p_value = enrichment_analysis( number_diff,number_all,len(all_spec_go), len(diff_sepc_go) )
	enrich_cache.append(each_go+'\t'+ go_def[  each_go ] +'\t%s\t%s\t%s'%( len(all_spec_go) ,len(diff_sepc_go),p_value  ))
	p_value_list.append(p_value)
q_value_list = iter(fdr(p_value_list))
ENRICHMENT = open( '%s.GO-mapping.Enrichment'%(output_result),'w'  )
ENRICHMENT.write("ROOT\tROOT_DEF\tTOTAL_Number\tDiff_Number\tP_value\tQ_value\n")
for key in enrich_cache:
	ENRICHMENT.write(key+'\t%s\n'%(q_value_list.next()))
MAPPING_DETAIL = open("%s.GO-mapping.detail"%(output_result),'w')
for each_id in id_go:
	MAPPING_DETAIL.write( each_id+'\t'+'\t'.join( id_go[ each_id ]  ) +'\n' )


