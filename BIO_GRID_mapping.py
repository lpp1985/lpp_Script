#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/26
from lpp import *
RAW = open( sys.argv[1],'rU' )
title = sys.argv[1].split('.')[0]
RAW.next()
END_ACC = open(title+'.acc' , 'w')
all_acc = {}
for line in RAW:
	line_l = line.split('\t')
	acc = re.search('gi\|(\d+)',line_l[6]).group(1)
	name = line_l[2].split()[0]
	END_ACC.write( name +'\t'+acc+'\n' )
ALL_links = open( '/home/wyh/Project/mydatabase/Biogrid/3.1.86/BIOGRID-IDENTIFIERS-3.1.86.tab.txt','rU' )
for line in ALL_links:
	if line.startswith('BIOGRID_ID'):
		break
gene2grid = {}
gene2taxon = {}
for line in ALL_links:
	line_l = line[:-1].split('\t')
	if 'GENBANK_PROTEIN_GI'== line_l[2]:
		gene2grid[ line_l[1]  ] = line_l[0]
		gene2taxon[ line_l[1]  ] = line_l[3]
RAW = open(END_ACC.name,'rU')
MAPPING = open( title+'.mapping'  ,'w')
MAPPING.write(  'ID\tGI\tGI-link\tBIOGRID_ID\tTAXON\n'  )
for line in RAW:
	line_l = line[:-1].split('\t')
	if line_l[-1] in gene2grid:
		MAPPING.write( line_l[0]+'\t'+line_l[1]+ '\t'+'http://www.ncbi.nlm.nih.gov/protein/%s'%( line_l[1] )+'\t'+gene2grid[  line_l[1]  ] +'\t'+gene2taxon[  line_l[1] ]+'\n' )

INTERACT = open( title+'.interact'  ,'w')
INTERACT.write(  'BIOGRID-ID1\tID\tBIOGRID-ID2\tID\n'   )

GRID = open( title+'.grid'  ,'w')
all_interact = {}
ALL_INT = open( '/home/wyh/Project/mydatabase/Biogrid/3.1.86/acc.out','rU'  )
all_need = {}
MAPPING_READ = open(MAPPING.name,'rU')
MAPPING_READ.next()
for line in MAPPING_READ:
	line_l = line[:-1].split('\t')
	all_need[  line_l[3] ] = line_l[0]
for line in ALL_INT:
	line_l = line.split('\t')
	cache=[]
	if line_l[0] in all_need or line_l[1] in all_need:
		for each_l in line_l[:2]:
			if each_l in all_need:
				cache.append( each_l+'\t'+  all_need[ each_l ] )
			else:
				cache.append( each_l+'\t-' )
		INTERACT.write( '\t'.join( cache )+'\n'  ) 
		all_interact[ '\t'.join( line_l[:2]  )  ]= ''
for each_key in all_interact:
	GRID.write( each_key+'\n' )
	
	