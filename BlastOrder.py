#!/usr/bin/python
#Author=LPP
from lpp import *
RAW = open(sys.argv[1],'rU')
location_record_hash = Ddict()
already = Ddict() 
def median(  List  ):
	location = sum( List ) /2
	return location
	
RAW.next()
already = {}
contig_order = {}
for line in RAW:
	line_l = line[:-1].split('\t')
	if line_l[2] in already:
		continue
	location_list = [ int(line_l[15]), int(line_l[16])    ]
	location =median( location_list)
	frame = line_l[18]
	contig_order[ location ] = line_l[2]+'\t'+frame
	already[ line_l[2] ] = ''
END = open( sys.argv[2],'w'  )
for key in sorted(  contig_order  ):
	END.write( '%s\t%s\n'%( key,contig_order[key]   )   )