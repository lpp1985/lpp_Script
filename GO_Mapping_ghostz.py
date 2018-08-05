#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/27
from lpp import *


RAW = open( sys.argv[1],'rU'  )

RAW.next()
geneid_id = Ddict()
for line in RAW:
	
	line_l = line.split('\t')

	name = line_l[0].split()[0]
	gi = re.search( '(?:sp|tr)\|(\w+)' ,line_l[1] ).group(1)

	output = '-\t-\t->\t'+name+'\t'+gi
	geneid_id[ gi ][ name ] = output

ALL_GO = open( '/home/lpp/Database/GO/gene_association.goa_uniprot','rU'  )
ALL_GO.next()
go_id = Ddict()
id_go = Ddict()
go_def = {}
DEF = open( '/home/lpp/Database/GO/NAME_DEF.list','rU'  )
change_go = File_Ddict( open(  '/home/lpp/Database/GO/relationship.alter','rU'  )  ).read(1,2)
for line in DEF:
	line_l = line.split('\t')
	go_def[ line_l[0] ] = line_l[1]

for line in ALL_GO:
	line_l = line.split('\t')
	if line_l[0] !='UniProtKB':
		continue
	if line_l[1] in geneid_id:

		for each_name in geneid_id[  line_l[1]  ]:
			#try:

				if line_l[4] not in change_go:
					
					go_id[ line_l[4] ][ each_name ] = geneid_id[  line_l[1]  ][ each_name  ]+'\t'+line_l[4]+'\t'+go_def[  line_l[4] ]
					id_go[ each_name ][ line_l[4] ] = ''
				else:
					for each_altered in change_go[ line_l[4] ]:
						go_id[ each_altered ][ each_name ] = geneid_id[  line_l[1]  ][ each_name  ]+'\t'+each_altered+'\t'+go_def[  each_altered ]
						id_go[ each_name ][ each_altered ] = ''
			#except:
				#pass
			
print(go_id)
son_father = File_Ddict( open('/home/lpp/Database/GO/relationship.son','rU' ) ).read(2,1)
all_stop = File_Ddict(  open('/home/lpp/Database/GO/ROOT.root' ,'rU') ).read(1)
root_static = Ddict()

def recall( leaf ):
	global root_static,need_check,path,go_id,number_check
	for each_son in  son_father[ leaf ]:
		if each_son not in all_stop:

			recall( each_son )
		else:
			root_static[ each_son ][ need_check ] = ''
			for each_id in go_id[need_check]:
				number_check[ each_son ][ each_id ] = ''
number_check = Ddict()
for each_go in go_id:
	need_check = each_go
	each_son = recall( each_go )

END = open( sys.argv[1]+'.GO-mapping.list','w'  )
END.write( 'ROOT\tROOT_DEF\tNumber\tID\tUnirpotKB\tGO\tGO_DEF\n'  )
for each_go in reversed( sorted(root_static,key= lambda x: len( number_check[ x  ]) ) ):
	END.write( each_go+'\t'+ go_def[  each_go ] +'\t%s\n'%( len( number_check[ each_go  ]  )  ) )
	for each_son in root_static[ each_go ]:
		for each_id in go_id[ each_son ]:
			
			END.write(  go_id[ each_son ][  each_id ]+'\n' )
			
			
			
			
END_DETAIL = open( sys.argv[1]+'.GO-mapping.detail','w'  )
for each_id in id_go:
	END_DETAIL.write( each_id+'\t'+'\t'.join( id_go[ each_id ]  ) +'\n' )
END_CAT = open( sys.argv[1]+'.GO-mapping.cat','w'  )	
for each_id in id_go:
	for key in id_go[ each_id ]:
        	END_CAT.write( each_id+'\t'+key +'\n' )		
			
			
