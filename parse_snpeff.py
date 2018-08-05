#!/usr/bin/python
#coding:utf-8
# Author:   --<>
# Purpose: 
# Created: 2013/5/18


from lpp import *
order_list  = []
has_order = {}
promoter = Ddict()
mutation = Ddict()
RAW = open( sys.argv[1],'rU'  )
END = open( sys.argv[2],'w' )
for line in RAW:
	if line.startswith("#"):
		continue
	
	try:
		all_eff = re.search( 'EFF=([^\t]+)',line  ).group(1)
	except:
		print(line)
		continue
	all_eff_list = all_eff.split(',')
	line_l = line.strip().split('\t')
	if 'INTERGENIC' in all_eff:
		all_regula_gene_list = filter(lambda x: x.startswith('UPSTREAM'), all_eff_list)
		effect_gene = []
		for each_gene in all_regula_gene_list:
			effect_gene.append(  each_gene.split('|')[8]  ) 
		if not  effect_gene:
			continue
		promoter_name = '_'.join( sorted( effect_gene) )+'-promoter'
		if promoter_name not in has_order:
			order_list.append( promoter_name )
		has_order[ promoter_name ] = ''
		promoter[ promoter_name ] [ '\t'.join( [  line_l[ 1 ], line_l[3],line_l[4]      ])   ] = ''
#	high_effect = []
#	for item in all_eff_list:
#		if 'MODIFIER' not in all_eff_list:
#			high_effect.append( item )
	high_effect = filter( lambda x: 'MODIFIER' not in x,all_eff_list   )
	for each_eff in high_effect:
		print(each_eff)
		category = re.search(  '(\S+)\(',each_eff ).group(1)
		all_info = re.search( '\((.+)\)',each_eff  ).group(1)
		all_info_list = all_info.split( '|' )
		level = all_info_list[0]
		nul_change = all_info_list[2]
		aa_change = all_info_list[3]
		gene_name = all_info_list[ 8 ]
		if gene_name not in has_order:
			order_list.append( gene_name )
		has_order[ gene_name  ] = ''
		mutation[  gene_name   ][   '\t'.join(  [  line_l[ 1 ], line_l[3],line_l[4]       , nul_change,aa_change,category,level ])    ] = ''
		
END.write(  '\t'.join(  ["GENE_ID","Genome Pos","Ref","Alt","Gene Seq Change","AA Change","Category","Level"  ]  )+'\n'  )

for each_gene in order_list:
	END.write(  each_gene )
	
		
	if each_gene.endswith('-promoter') and each_gene in promoter:

		for each_key2 in promoter[each_gene]:
			
			END.write(  '\t'+each_key2+'\t-'*4+'\n'  )
			
	else:
		for each_key2 in mutation[ each_gene   ]:
			END.write( '\t'+each_key2+'\n'  )
