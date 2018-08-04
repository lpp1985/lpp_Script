#!/usr/bin/env python
#coding:utf-8
# Author:   --<>
# Purpose: 
# Created: 2013/6/8
# Convert maq's result to vcf
# Usage $SCRIPT snp indelpe OUTPUT
from lpp import *
SNP = open( sys.argv[1],'rU'   )
INDEL = open( sys.argv[2],'rU'   )
END = open( sys.argv[3],'w'  )
END.write( '''
##fileformat=VCFv4.1
##fileDate=20130606
##reference=NC_000913
##INFO=DP,1,Integer,"Total Depth of Coverage"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
'''
           )

for line in SNP:
	line_l = line.split( '\t' )
	new_cache = [line_l[0],line_l[1],'.',line_l[2],line_l[3],line_l[4],'PASS','DP=%s'%( line_l[5] )]
	END.write( '\t'.join(new_cache )+'\n')
for line in INDEL:
	line_l = line.split( '\t'  )
	seq_cache = re.search( '\:(\w+)',line_l[4]  ).group(1)
	if line_l[4].startswith('-'):
		old_seq = line_l[4][-1]+seq_cache
		new_seq = line_l[4][-1]
	else:
		new_seq = line_l[4][-1]+seq_cache
		old_seq = line_l[4][-1]
	new_cache = [ line_l[0], '%s'%( int(line_l[1])-1 ),'.',  old_seq,new_seq ,'10','PASS','DP=%s'%( int( line_l[3] ) - int( line_l[9] )  )   ]
	END.write( '\t'.join( new_cache )+'\n' )
	
	
	
	
