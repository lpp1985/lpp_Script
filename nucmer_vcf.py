#!/usr/bin/python
#coding:utf-8
# Author:   --<>
# Purpose: 
# Created: 2014/4/24

from lpp import *
REF = fasta_check(  open( sys.argv[1],'rU' ) )
ref_seq = {}
for t,s in REF:
	ref_seq[ t.strip()[1:].split()[0]  ] = re.sub( "\s+",'',s )
	
SNP = open( sys.argv[2],'rU' )
all_change = Ddict()


lastposion = 0

for line in SNP:
	if '[' in line:
		continue
	line_l = line.split("\t")
	position,ref,change,refID = int(line_l[0])-1,line_l[1],line_l[2],line_l[-2]
	if change=='.':
		change = ''
	if ref =='.':
		ref=''

	if position==lastposion+1:
		position = lastposion
			
			
			
	if refID in all_change and position in all_change[ refID ] and "que" in all_change[ refID ] [position] :
		
		
		
		all_change[ refID ] [position] ["que"] += change
	else:
		all_change[ refID ] [position] ["que"] = change
		
	if refID in all_change and position in all_change[ refID ] and "ref" in all_change[ refID ] [position] :
		
		all_change[ refID ] [position] ["ref"] += ref
	else:
		all_change[ refID ] [position] ["ref"] = ref	
	
	
	lastposion = int(line_l[0])-1
		

END = open( sys.argv[3],'w' )
END.write( """##fileformat=VCFv4.1
##fileDate=20130702
##reference=NC_000913
##INFO=DP,1,Integer,"Total Depth of Coverage"
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
"""  )
for ref in all_change:
	for key2 in sorted( all_change[ref] ):
		raw_base = ref_seq[ ref  ][ key2-1  ]
		END.write(  
		    "%s\t%s\t.\t%s\t%s"%(
		        ref,
		        key2,
		        raw_base+ all_change[ ref  ][ key2  ][ "ref"  ] ,
		        raw_base+all_change[ ref  ][ key2  ][ "que"  ]
		    )
		    +"\t100.00\tPASS\tDP=100\n"
		)
